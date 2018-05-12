#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <jsoncpp/json/json.h>
#include <map>
#include <vector>

#pragma once

namespace structan
{

struct joint_t
{
    double x, y, z;
};

struct element_t
{
    int j1, j2;
};

struct eqno_t
{
    int u1, u2, u3, r1, r2, r3;
};

class solver_t
{
    std::vector<joint_t> joints;
    std::vector<element_t> elements;
    std::map<int, eqno_t> eqnos; // a map from joint no/dof to eq no
    std::vector<double> k, m, u, u_, u__, k_, r;
    double dt, alpha, delta, a0, a2, a3, a6, a7;
    bool fx, fy, fz;
    std::map<double, double> ax, ay, az;
    double time;

    void read_joints(std::string fname)
    {
        std::ifstream file(fname);
        std::string dummy;

        // ignore first three lines
        std::getline(file, dummy);
        std::getline(file, dummy);
        std::getline(file, dummy);

        while (true)
        {
            joint_t j;
            file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> j.x >> j.y >> j.z;
            if (file.eof())
                break;
            joints.push_back(j);
        }

        file.close();
    }

    void read_elements(std::string fname)
    {
        std::ifstream file(fname);
        std::string dummy;

        // ignore first three lines
        std::getline(file, dummy);
        std::getline(file, dummy);
        std::getline(file, dummy);

        element_t e;
        while (true)
        {
            file >> dummy >> e.j1 >> e.j2 >> dummy >> dummy >> dummy >> dummy >> dummy;
            if (file.eof())
                break;

            --e.j1;
            --e.j2;

            elements.push_back(e);
        }

        file.close();
    }

    void read_eqnos(std::string fname)
    {
        std::ifstream file(fname);

        // ignore first line
        std::string dummy;
        std::getline(file, dummy);

        while (true)
        {
            int n;
            eqno_t q;
            file >> n >> q.u1 >> q.u2 >> q.u3 >> q.r1 >> q.r2 >> q.r3;
            if (file.eof())
                break;
            --n;
            --q.u1;
            --q.u2;
            --q.u3;
            --q.r1;
            --q.r2;
            --q.r3;
            eqnos[n] = q;
        }

        file.close();
    }

    void read_k(std::string fname)
    {
        int n = m.size();
        k.resize(n * n);

        std::ifstream file(fname);

        // ignore first line
        std::string dummy;
        std::getline(file, dummy);

        int r, c;
        double v;
        while (true)
        {
            file >> r >> c >> v;
            if (file.eof())
                break;

            --r;
            --c;

            k[r * n + c] = v;
        }

        file.close();
    }

    void read_m(std::string fname)
    {
        std::string dummy;
        std::ifstream file(fname);

        // ignore first line
        std::getline(file, dummy);

        double v;
        while (true)
        {
            file >> dummy >> dummy >> v;
            if (file.eof())
                break;

            m.push_back(v);
        }

        file.close();
    }

    static std::map<double, double> read_timeseries(const std::string fname)
    {
        std::ifstream file(fname);
        std::map<double, double> result;
        double t, v;
        while (true)
        {
            file >> t >> v;
            if (file.eof())
                break;
            result[t] = v;
        }
        file.close();
        return result;
    }

    static double interpolate(std::map<double, double> &data, double t)
    {
        auto it = data.upper_bound(t);
        if (it == data.end())
            return (--it)->second;
        if (it == data.begin())
            return it->second;

        auto l = it;
        --l;

        const double delta = (t - l->first) / (it->first - l->first);
        return delta * it->second + (1 - delta) * l->second;
    }

    // http://www.sci.utah.edu/~wallstedt/LU.htm
    // Cholesky requires the matrix to be symmetric positive-definite
    static void Cholesky(int d, std::vector<double> &S, std::vector<double> &D)
    {
        for (int k = 0; k < d; ++k)
        {
            double sum = 0.;
            for (int p = 0; p < k; ++p)
                sum += D[k * d + p] * D[k * d + p];
            D[k * d + k] = std::sqrt(S[k * d + k] - sum);
            for (int i = k + 1; i < d; ++i)
            {
                double sum = 0.;
                for (int p = 0; p < k; ++p)
                    sum += D[i * d + p] * D[k * d + p];
                D[i * d + k] = (S[i * d + k] - sum) / D[k * d + k];
            }
        }
    }

    void solveCholesky(int d, std::vector<double> &LU, std::vector<double> &b, std::vector<double> &x)
    {
        double y[d];
        for (int i = 0; i < d; ++i)
        {
            double sum = 0.;
            for (int k = 0; k < i; ++k)
                sum += LU[i * d + k] * y[k];
            y[i] = (b[i] - sum) / LU[i * d + i];
        }
        for (int i = d - 1; i >= 0; --i)
        {
            double sum = 0.;
            for (int k = i + 1; k < d; ++k)
                sum += LU[k * d + i] * x[k];
            x[i] = (y[i] - sum) / LU[i * d + i];
        }
    }

  public:
    solver_t(std::string fname)
    {
        Json::Reader reader;
        Json::Value root;
        std::ifstream file(fname);
        reader.parse(file, root, false);
        file.close();

        read_joints(root["joints_file"].asString());
        read_elements(root["elements_file"].asString());
        read_eqnos(root["eqnos_file"].asString());
        read_m(root["m_file"].asString());
        read_k(root["k_file"].asString());

        alpha = root["alpha"].asDouble();
        delta = root["delta"].asDouble();

        fx = root["fx"].asBool();
        fy = root["fy"].asBool();
        fz = root["fz"].asBool();

        if (fx)
            ax = read_timeseries(root["ax_file"].asString());
        else
            ax[0] = 0;
        if (fy)
            ay = read_timeseries(root["ay_file"].asString());
        else
            ay[0] = 0;
        if (fz)
            az = read_timeseries(root["az_file"].asString());
        else
            az[0] = 0;

        u.resize(m.size());
        u_.resize(m.size());
        u__.resize(m.size());

        r.resize(m.size());
        k_.resize(k.size());

        a0 = 1 / (alpha * dt * dt);
        a2 = 1 / (alpha * dt);
        a3 = 1 / (2 * alpha) - 1;
        a6 = dt * (1 - delta);
        a7 = delta * dt;

        int neq = m.size();
        std::copy_n(k.begin(), neq, k_.begin());
        for (int i = 0; i < neq; ++i)
            k_[i * neq + i] += a0 * m[i];

        Cholesky(neq, k_, k_);

        time = 0;
    }

    void step()
    {
        time += dt;
        int neq = m.size();
        std::vector<double> uu(neq), uu_(neq), uu__(neq);

        std::fill_n(r.begin(), neq, 0);
        for (const auto &i : eqnos)
        {
            r[i.second.u1] = interpolate(ax, time) * m[i.second.u1];
            r[i.second.u2] = interpolate(ay, time) * m[i.second.u2];
            r[i.second.u3] = interpolate(az, time) * m[i.second.u3];
        }
        for (int i = 0; i < neq; ++i)
            r[i] += m[i] * (a0 * u[i] + a2 * u_[i] + a3 * u__[i]);

        solveCholesky(neq, k_, r, uu);

        for (int i = 0; i < neq; ++i)
        {
            uu__[i] = a0 * (uu[i] - u[i]) - a2 * u_[i] - a3 * u__[i];
            u_[i] = u_[i] + a6 * u__[i] + a7 * uu__[i];
        }

        u = uu;
        u_ = uu_;
        u__ = uu__;
    }

    void write(std::string fname)
    {
        int nnodes = joints.size();
        int nelems = elements.size();

        std::ofstream file(fname);

        file << "# vtk DataFile Version 2.0" << std::endl
             << "structure data" << std::endl
             << "ASCII" << std::endl
             << "DATASET UNSTRUCTURED_GRID" << std::endl;

        file << "POINTS " << nnodes << " double" << std::endl;
        for (int i = 0; i < nnodes; ++i)
        {
            auto d = eqnos[i];
            file << joints[i].x + u[d.u1] << " "
                 << joints[i].y + u[d.u2] << " "
                 << joints[i].z + u[d.u3] << std::endl;
        }

        file << "CELLS " << nelems << " " << 3 * nelems << std::endl;
        for (int i = 0; i < nelems; ++i)
            file << "2 " << elements[i].j1 << " " << elements[i].j2 << std::endl;

        file << "CELL_TYPES " << nelems << std::endl;
        for (int i = 0; i < nelems; ++i)
            file << 3 << std::endl;

        file.close();
    }
};
}
