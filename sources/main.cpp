#include <cmath>
#include <functional>
#include <iostream>

#include <eigen3/Eigen/Dense>

int main(int argc, char const *argv[]) {
    // test potential
    auto potential = [](double r) {
        const double a[20] = {-0.065, 15939.056, -29646.778000000002, -6269.777, 44952.358,
                            8709.016, -100549.29, 594784.152, -995239.126, -1.14496717e7,
                            4.606463055e7, 3.7466657300000004e7, -5.439157146e8,
                            9.36483394e8, 1.387879748e9, -8.400905473e9, 1.5781752106e10,
                            -1.5721037673e10, 8.376043061000001e9, -1.8898488e9};

        const double bt     = -0.17;
        const double Rm     = 4.6719018;
        const double Tm     = -1081.6384;
        const double Ri     = 3.963;
        const double nnn    = 12.362;
        const double AAA    = -1.3328825e3;
        const double BBB    = 3.321662099e10;
        const double Ra     = 10.5;
        const double C62    = 1.525e7;
        const double C8     = 5.159e8;
        const double C10    = 1.91e10;
        const double tobohr = 0.52917721092;
        const double tocm   = 219474.6314;

        r *= tobohr;

        using std::pow;

        if (r < Ri)
            return (AAA + BBB / pow(r, nnn)) / tocm;
        else if (r < Ra) {
            double res = Tm;
            for (int i = 0; i < 20; ++i)
                res += a[i] * pow((r - Rm) / (r + bt * Rm), i + 1);
            return res / tocm;
        } else
            return (-C62 / pow(r, 6) - C8 / pow(r, 8) - C10 / pow(r, 10)) / tocm;
    };

    const int npoints = 4000;
    const double rmin = 4.0;
    const double rmax = 10;

    const double m = 80121.07031084;
    const int j    = 0;

    auto centrifugal = [m, j](double r) {
        return j * (j + 1) / 2. / m / pow(r, 2);
    };

    const double npo  = npoints + 1;
    const double span = rmax - rmin;
    const double dx   = span / npo;

    using Eigen::MatrixXd;

    MatrixXd T(npoints, npoints);
    MatrixXd V(npoints, npoints);
    T.setZero();
    V.setZero();

    using std::sin;

    for (int i = 0; i < npoints; ++i)
        for (int j = 0; j <= i; ++j) {
            if (i == j) {
                V(i, j) = potential(rmin + (i+1) * dx) + centrifugal(rmin + (i+1) * dx);
                T(i, j) = pow(M_PI, 2) / 4. / pow(span, 2) * ((2. * pow(npo, 2) + 1) / 3. - 1. / pow(sin(M_PI * (i+1) / npo), 2));
            }
            else {
                T(i, j) = pow(-1, i - j) * pow(M_PI, 2) / 4. / pow(span, 2) * 
                ( 1. / pow(sin(M_PI * (i - j) / 2. / npo ), 2) -  1. / pow(sin(M_PI * (i + j + 2) / 2. / npo ), 2) );
                T(j, i) = T(i, j);
            }

        }

    MatrixXd H = T + V;

    Eigen::SelfAdjointEigenSolver<MatrixXd> es(H);

    const double tomhz = 6.579695e9;

    using Eigen::VectorXd;

    VectorXd lowest = es.eigenvalues().head(4);
    std::cout << '\n';
    std::cout << V.diagonal() << '\n';
    std::cout << '\n';
    std::cout << lowest << '\n';
    std::cout << '\n';

}
