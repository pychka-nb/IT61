#include <iostream>
//#include <time.h>
//#include <math.h>
#include <fstream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
void svd_hestenes(int m_m, int n_n, float* a, float* u, float* v, float* sigma)
{
	float thr = 0.000001f;
	int n, m, i, j, l, k, lort, iter, in, ll, kk;
	float alfa, betta, hamma, eta, t, cos0, sin0, buf, s, r;
	n = n_n;
	m = m_m;
	for (i = 0; i < n; i++)
	{
		in = i * n;
		for (j = 0; j < n; j++)
			if (i == j) v[in + j] = 1.;
			else v[in + j] = 0.;
	}
	for (i = 0; i < m; i++)
	{
		in = i * n;
		for (j = 0; j < n; j++)
		{
			u[in + j] = a[in + j];
		}
	}

	iter = 0;
	while (1)
	{
		lort = 0;
		iter++;
		for (l = 0; l < n - 1; l++)
			for (k = l + 1; k < n; k++)
			{
				alfa = 0.; betta = 0.; hamma = 0.;
				for (i = 0; i < m; i++)
				{
					in = i * n;
					ll = in + l;
					kk = in + k;
					alfa += u[ll] * u[ll];
					betta += u[kk] * u[kk];
					hamma += u[ll] * u[kk];
				}

				if (sqrt(alfa * betta) < 1.e-10)	continue;
				if (fabs(hamma) / sqrt(alfa * betta) < thr)
					continue;

				lort = 1;
				eta = (betta - alfa) / (2.f * hamma);
				t = (eta / fabs(eta)) /
					(fabs(eta) + (float)sqrt(1.f + eta * eta));
				cos0 = 1.f / (float)sqrt(1.f + t * t);
				sin0 = t * cos0;

				for (i = 0; i < m; i++)
				{
					in = i * n;
					buf = u[in + l] * cos0 - u[in + k] * sin0;
					u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
					u[in + l] = buf;

					if (i >= n) continue;
					buf = v[in + l] * cos0 - v[in + k] * sin0;
					v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
					v[in + l] = buf;
				}
			}

		if (!lort) break;
	}

	for (i = 0; i < n; i++)
	{
		s = 0.;
		for (j = 0; j < m; j++)	s += u[j * n + i] * u[j * n + i];
		s = (float)sqrt(s);
		sigma[i] = s;
		if (s < 1.e-10)	continue;
		for (j = 0; j < m; j++)	u[j * n + i] = u[j * n + i] / s;
	}
	for (i = 0; i < n - 1; i++)
		for (j = i; j < n; j++)
			if (sigma[i] < sigma[j])
			{
				r = sigma[i]; sigma[i] = sigma[j]; sigma[j] = r;
				for (k = 0; k < m; k++)
				{
					r = u[i + k * n]; u[i + k * n] = u[j + k * n]; u[j + k * n] = r;
				}
				for (k = 0; k < n; k++)
				{
					r = v[i + k * n]; v[i + k * n] = v[j + k * n]; v[j + k * n] = r;
				}
			}

	//return iter;
}

int main()
{
    srand(time(NULL));
    /*int n, m;

    cout << "Vvedite razmer m,n" << endl;
    cin >> m >>n;
	float** a;
	float* x = new float[m];
	float* xNew = new float[m];
	float* b = new float[m];
	float* a1 = new float[m * n];
	float* aTr = new float[m * n];
	float* aNew = new float[m * n];
	float* u = new float[m * m];
	float* vTr = new float[n * n];
	float* uTr = new float[m * m];
	float* v = new float[n * n];
	float* sigma = new float[m * n];
    a = new float* [m];
    for (int i = 0; i < m; i++) a[i] = new float[n];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = rand() / (1.0*RAND_MAX);
            x[j]= rand() / (1.0 * RAND_MAX);
            b[j] = 0.;
			xNew[j] = 0.;
			aTr[i * m + j] = 0.;
			aNew[i * m + j] = 0.;
            //cout << a[i][j]<<" ";
        }
        //cout << endl;
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            b[i] += a[i][j] * x[j];
			a1[m * i + j] = a[i][j];
        }
    }
	svd_hestenes(m, n, a1, u, vTr, sigma);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			uTr[i + j * m] = u[i * m + j];
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			v[i + j * n] = vTr[i * n + j];
		}
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				sigma[i * m + j] = 1.0 / sigma[i * m + j];
			}
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			for (int k = 0; k < n; k++) {
				aTr[i * m + j] += v[i * n + k] * sigma[k * m + j];
			}
		}
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			for (int k = 0; k < m; k++) {
				aNew[i * m + j] += aTr[i * m + k] * uTr[k * m + j];
			}
		}
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			xNew[i] += aNew[m*i+j] * b[j];
			
		}
	}
	for (int i = 0; i < m; i++) {
		cout << x[i] << " " << xNew[i] << endl;
	}
	/*for (int i = 0; i < n * n; i++) {
		cout << vTr[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < n * n; i++) {
		cout << v[i] << " ";
	}*/

    /*for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
    for (int i = 0; i < m; i++) {
        cout << x[i] << " " << b[i] << endl;
    }*/
	int n, m;

	cout << "Vvedite razmer m,n" << endl;
	cin >> m >> n;
	MatrixXd A = MatrixXd::Random(m, n).array().abs();
	VectorXd X = VectorXd::Random(n).array().abs();
	VectorXd B = A * X;
	JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
	MatrixXd U = svd.matrixU();
	MatrixXd V = svd.matrixV();
	VectorXd singVal = svd.singularValues();
	MatrixXd singValInv = singVal.asDiagonal().inverse();
	MatrixXd Aps = V * singValInv * U.transpose();
	VectorXd Xps = Aps * B;
	double err1 = (X - Xps).norm();
	double err2 = (A * Xps - B).norm();

	ofstream outputFile("output.txt");
	if (outputFile.is_open()) {
		outputFile << "Матрица A:\n" << A << endl << endl;
		outputFile << "Матрица A#:\n" << Aps << endl << endl;
		outputFile << "Сингулярные значения:\n" << singVal << endl << endl;
		outputFile << "Вектор X:\n" << X << endl << endl;
		outputFile << "Вычисленный вектор X#:\n" << Xps << endl << endl;
		outputFile << "Невязка |X - X#|: " << err1 << endl;
		outputFile << "Невязка |A * X# - B|: " << err2 << endl;
		outputFile.close();
		cout << "Finish" << endl;
	}
	else {
		cerr << "Ошибка открытия файла output.txt" << endl;
		return 1;
	}

	return 0;
}



