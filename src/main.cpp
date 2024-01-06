#include "Eigen/Dense"
#include <iostream>
#include <omp.h>

#include <viewer.hpp>

using namespace std;

int main() {
    Viewer viewer;
    std::vector<Eigen::Vector2d> points;
    const bool &is_closed = viewer.is_closed();

    auto reset_buffer = [&viewer](const std::vector<Eigen::Vector2d> &points) {
        viewer.reset_buffer(points);
    };
    auto update_buffer = [&viewer](const std::vector<Eigen::Vector2d> &points) {
        viewer.update_buffer(points);
    };

#pragma omp parallel sections num_threads(2) default(shared)
    {
#pragma omp section
        { viewer.show(); }
#pragma omp section
        {
            double t = 0;
            Eigen::Vector2d x(sin(t) / 2 + 0.5, cos(t) / 2 + 0.5);
            std::vector<Eigen::Vector2d> nodes_tmp;
            nodes_tmp.push_back(x);
            reset_buffer(nodes_tmp);

            while (!is_closed) {
                x = Eigen::Vector2d(sin(t) / 2 + 0.5, cos(t) / 2 + 0.5);
                t += 1e-7;

                nodes_tmp.clear();
                nodes_tmp.push_back(x);
                update_buffer(nodes_tmp);
            }
        }
    }
    return 0;
}