using System;
using System.Linq;

namespace QR_iter
{
    class Program
    {
        static void Main(string[] args)
        {
            var slau = SLAU.GetSLAU();
            slau.PrintSLAU();
            Console.WriteLine("Выберите метод:\n1 - решение с помощью QR разложения;\n2 - метод простой итерации;\n3 - LU разложение;\n" +
                              "4 - метод квадратных корней;\n5 - нахождение минимального собственного числа степенным методом;\n" +
                              "6 - нахождение собственных значений QR-алгоритмом");
            switch (int.Parse(Console.ReadLine()))
            {
                case 1:
                    Console.WriteLine();
                    slau.SolveSlauQR(slau.Matrix, slau.B); //v
                    break;
                case 2:
                    Console.WriteLine();
                    slau.IterativeMethod();  //хуйня
                    break;
                case 3:
                    Console.WriteLine();
                    slau.LUdecomposition(); //v
                    break;
                case 4:
                    Console.WriteLine();
                    slau.UUTdecomposition(); //v
                    break;
                case 5:
                    Console.WriteLine();
                    slau.PowerIteration();  //v
                    break;
                case 6:
                    Console.WriteLine();
                    slau.QRdecomposition(); //числа
                    break;
                default:
                    Console.WriteLine();
                    Console.WriteLine("Неверный номер метода");
                    break;
            }


            Console.ReadLine();
        }
    }

    public class SLAU //экземпляр класса — СЛАУ (все коэффициенты и свободные члены)
    {
        
        //конструктор принимает количество уравнений
        public SLAU(int n)
        {
            Matrix = new double[n, n];
            B = new double[n];
            N = n;
        }

        //ввод данных с клавиатуры
        public static SLAU GetSLAU()  
        {
            Console.WriteLine("Введите количество уравнений:");
            var n = int.Parse(Console.ReadLine());
            var slau = new SLAU(n);
            Console.WriteLine();
            Console.WriteLine("Введите коэффициенты (через пробел, разделяя уравнения переносом строки):");
            for (int i = 0; i < n; i++)
            {
                var str = Console.ReadLine().Split(' ');
                for (int j = 0; j < n; j++)
                {
                    slau.Matrix[i, j] = double.Parse(str[j]);
                }
            }
            Console.WriteLine();
            Console.WriteLine("Введите свободные члены (через пробел):");
            var str2 = Console.ReadLine().Split(' ');
            for (int j = 0; j < n; j++)
            {
                slau.B[j] = double.Parse(str2[j]);
            }
            Console.WriteLine();
            return slau;
        }

        //вывод коэффициентов и свободных членов
        public void PrintSLAU() 
        {
            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(0); j++)
                {
                    Console.Write("{0,5}", Matrix[i, j]);                 
                }
                Console.Write(" | {0,5}", B[i]);
                Console.WriteLine();
                Console.WriteLine();
            }
        }

        //вывод решения СЛАУ
        public void PrintResult(double[] res)
        {
            for (int i = 0; i < res.Length; i++)
            {
                Console.WriteLine("{0}) {1:0.000}", i+1, res[i]);
            }
            Console.WriteLine();
            Console.WriteLine();
        }

        //вывод матрицы
        public void PrintMatrix(double[,] matrix)
        {
            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(0); j++)
                {
                    Console.Write("{0,7:0.000}  ", matrix[i, j]);
                }
                Console.WriteLine();
                Console.WriteLine();
            }
        }




        //МЕТОД ПРОСТЫХ ИТЕРАЦИЙ
        //реализация метода простых итераций
        public void IterativeMethod()
        {
            Console.WriteLine("Задайте точность вычислений для итерационного метода:");
            var eps = double.Parse(Console.ReadLine());
            Console.WriteLine();
            var c = new double[N];
            var xk = B; //x(k)
            var xknext = new double[N]; // x(k+1)
            var matrixB = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i != j)
                    {
                        matrixB[i, j] = -Matrix[i, j] / Matrix[i, i];
                    }
                    else
                    {
                        matrixB[i, j] = 0;
                    }

                }
                c[i] = B[i] / Matrix[i, i];
            }
            if (!CheckIfMatrixIsValid(matrixB))
            {
                Console.WriteLine("Не выполнено достаточное условие сходимости итерационной последовательности");
            }
            else
            {
                while (!CheckConditionForIterativeMethod(xk, xknext, FindMatrixNorm(matrixB), eps))
                {
                    for (int i = 0; i < N; i++)
                    {
                        xknext[i] = 0;
                        for (int j = 0; j < N; j++)
                        {
                            xknext[i] += matrixB[i, j] * xk[j];
                        }
                        xknext[i] += c[i];
                    }
                    xk = xknext;
                }
                Console.WriteLine("Решение системы:");
                PrintResult(xknext);
            }
        }

        //расчет бесконечной нормы матрицы
        private double FindMatrixNorm(double[,] matrix)
        {
            var sum = new double[N];
            for (int i = 0; i < N; i++)
            {
                sum[i] = 0;
                for (int j = 0; j < N; j++)
                {
                    sum[i] += Math.Abs(matrix[i, j]);
                }
            }
            return sum.Max();
        }

        //проверка достаточного условия сходимости итерационной последовательности
        private bool CheckIfMatrixIsValid(double[,] matrix)
        {
            double Ek = 0;
            for (int i = 0; i < B.Length; i++)
            {
                for (int j = 0; j < B.Length; j++)
                {
                    Ek = Ek + matrix[i, j] * matrix[i, j];
                }
            }
            Ek = Math.Sqrt(Ek);
            if (Ek > 1)
            {
                return false;
            }
            return true;
        }

        //проверка выполнения условия для окончания повторения итераций (метод простых итераций)
        private bool CheckConditionForIterativeMethod(double[] xk, double[] xknext, double q, double eps)
        {
            var x = new double[N];
            for(int i = 0; i < N; i++)
            {
                x[i] = Math.Abs(xknext[i] - xk[i]);
            }
            if (x.Max() <= (1 - q) * eps / q ) { return true;}
            return false;
        }




        //МЕТОД LU-РАЗЛОЖЕНИЯ
        //реализация метода LU-разложений
        public void LUdecomposition()
        {
            var n = Matrix.GetLength(0);
            var L = new double[n, n];
            var U = new double[n, n];

            GetLU(L, U, n, Matrix);
            
            var trMatrix = new double[n, n];
            var sum1 = 0.0;
            for (int i = n-1; i >= 0; i--)
            {
                for (int j = n-1; j >= 0; j--)
                {
                    if (i == j)
                    {
                        for (int k = j + 1; k < n; k++)
                        {
                            sum1 += U[j, k] * trMatrix[k, j];
                        }
                        trMatrix[j, j] = (1 - sum1) / U[j, j];
                        sum1 = 0;
                    }
                    else if (i < j)
                    {
                        for (int k = i + 1; k < n; k++)
                        {
                            sum1 += U[i, k] * trMatrix[k, j];
                        }
                        trMatrix[i, j] = -sum1 / U[i, i];
                        sum1 = 0;
                    }
                    else
                    {
                        for (int k = j + 1; k < n; k++)
                        {
                            sum1 += U[k, j] * trMatrix[i, k];
                        }
                        trMatrix[i, j] = -sum1;
                        sum1 = 0;
                    }
                }
            }
            var invMatrix = GetInvMatrix();
            Console.WriteLine("Обратная матрица:");
            PrintMatrix(invMatrix);

            var x = SolveSlau(L, U, B);
            Console.WriteLine("Решение системы:");
            PrintResult(x);
        }

        //получение матриц L и U
        private void GetLU(double[,] L, double[,] U, int n, double[,] matrix)
        {
            for (int i = 0; i < n; i++)
            {
                U[0, i] = matrix[0, i];
                L[i, 0] = matrix[i, 0] / U[0, 0];
            }
            for (int i = 1; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    double sum = 0;
                    for (int k = 0; k < i; k++)
                    {
                        sum += L[i, k] * U[k, j];
                    }
                    U[i, j] = matrix[i, j] - sum;
                    {
                        sum = 0;
                        for (int k = 0; k < i; k++)
                        {
                            sum += L[j, k] * U[k, i];
                        }
                        L[j, i] = (matrix[j, i] - sum) / U[i, i];
                    }
                }
            }
        }

        private double[] SolveSlau(double[,] Matrix1, double[,] Matrix2, double[] B)
        {
            double[] y = new double[N];
            double[] x = new double[N];

            for (int i = 0; i < N; i++)   //прямая подстановка
            {
                double sum = 0;
                for (int j = 0; j < N; j++)
                {
                    sum += Matrix1[i, j] * y[j];
                }
                y[i] = (B[i] - sum) / Matrix1[i, i];
            }

            for (int i = N - 1; i >= 0; i--)   //Обратная подстановка
            {
                double sum = 0;
                for (int j = N - 1; j > i; j--)
                {
                    sum += Matrix2[i, j] * x[j];
                }
                x[i] = (y[i] - sum) / Matrix2[i, i];
            }
            return x;
        }




        //МЕТОД QR-РАЗЛОЖЕНИЯ
        //QR-алгоритм
        public void QRdecomposition()
        {
            var A = new double[B.Length, B.Length];
            var newA = new double[B.Length, B.Length];
            var Q = new double[B.Length, B.Length];
            var R = new double[B.Length, B.Length];
            var eigenvectors = new double[N,N];
            var eigenvalues = new double[N];
            var dif = 0.0;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    A[i, j] = Matrix[i, j];
                    if (i == j)
                    {
                        eigenvectors[i, j] = 1.0;
                    }
                    else
                    {
                        eigenvectors[i, j] = 0.0;
                    }
                }
            }

            do
            {

                GetQR(A, Q, R);

                var ev = MultiplyMatrix(eigenvectors, Q);
                eigenvectors = ev;

                newA = MultiplyMatrix(R, Q);


                dif = 0.0;
                for (int i = 0; i < N; i++)
                {
                    if (dif < Math.Abs(eigenvalues[i] - newA[i, i]))
                    {
                        dif = Math.Abs(eigenvalues[i] - newA[i, i]);
                    }
                    eigenvalues[i] = newA[i, i];
                }
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        A[i, j] = newA[i, j];
                    }
                }
            } while (dif > 0.000001);

            Console.WriteLine("Собственные числа:");
            PrintResult(eigenvalues);
            Console.WriteLine("Собственные векторы:");
            PrintMatrix(Rationing(eigenvectors));
        }

        //нормирование
        private double[,] Rationing(double[,] eigenvectors)
        {
            for (int j = 0; j < N; j++)
            {
                var sum = 0.0;
                for (int i = 0; i < N; i++)
                {
                    sum += eigenvectors[i, j] *eigenvectors[i, j];
                }
                for (int i = 0; i < N; i++)
                {
                    eigenvectors[i, j] = eigenvectors[i, j] / Math.Sqrt(sum);
                }
            }
            return eigenvectors;
        }

        //решение СЛАУ с помощью QR-разложения
        public double[] SolveSlauQR(double[,] matr, double[] b)
        {
            var matrixB = matr;
            var Q = new double[B.Length, B.Length];
            var R = new double[B.Length, B.Length];

            GetQR(matrixB, Q, R);

            var y = MultiplyMatrixVector(Q, b);
            var x = new double[N];

            x[N - 1] = y[N - 1] / R[N - 1, N - 1];
            for (int i = N - 2; i >= 0; i--)
            {
                x[i] = y[i];
                for (int j = N - 1; j > i; j--)
                {
                    x[i] -= R[i, j] * x[j];
                }
                x[i] /= R[i, i];

            }
            if (matr != Matrix)
            {
                return x;
            }
            else
            {
                Console.WriteLine("Решение системы:");
                PrintResult(x);
                return x;
            }
        }

        //получение матриц Q и R
        private void GetQR(double[,] A, double[,] Q, double [,] R)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    Q[i, j] = 0.0;
                    R[i, j] = 0.0;
                }
            }
            var QT = new double[N, N];
            for (int i = 0; i < N; i++)
                QT[i, 0] = A[i, 0];

            R[0, 0] = 1;
            var sum = 0.0;
            for (int k = 0; k < N; k++)
            {                
                for (int i = 0; i < N; i++)
                {
                    sum += A[i, k] * A[i, k];
                }
                R[k, k] = Math.Sqrt(sum);
                sum = 0;
                for (int i = 0; i < N; i++)
                    QT[i, k] = A[i, k] / R[k, k];

                for (int j = k; j < N; j++)
                {
                    for (int i = 0; i < N; i++)
                    {
                        sum += QT[i, k] * A[i, j];
                    }
                    R[k, j] = sum;
                    sum = 0;
                    for (int i = 0; i < N; i++)
                        A[i, j] = A[i, j] - QT[i, k] * R[k, j];
                }
            }
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    Q[i, j] = QT[i, j];
                }
            }
        }

        //перемножение матриц
        private double[,] MultiplyMatrix(double[,] R, double[,] Q)
        {
            var A = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    A[i, j] = 0;
                    for (int k = 0; k < N; k++)
                    {
                        A[i, j] += R[i, k] * Q[k, j];
                    }
                }
            }
            return A;
        }

        //умножение матрицы на вектор
        private double[] MultiplyMatrixVector(double[,] R, double[] B)
        {
            var A = new double[N];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    var sum = 0.0;
                    for (int k = 0; k < N; k++)
                    {
                        sum += R[i, k] * B[k];
                    }
                    A[i] = sum;
                }
            }
            return A;
        }



        //МЕТОД КВАДРАТНЫХ КОРНЕЙ
        //разложение на U и U^T
        public void UUTdecomposition()
        {
            var U = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                for (int j = 1; j < N; j++)
                {
                    var sum = 0.0;
                    for (int k = 0; k < i; k++)
                    {
                        sum += U[k, i] * U[k, i];
                    }
                    U[i, i] = Math.Sqrt(Matrix[i, i] - sum);
                    if (i != j)
                    {
                        sum = 0;
                        for (int k = 0; k < i; k++)
                        {
                            sum += U[k, i] * U[k, j];
                        }
                        if (j < i) U[i, j] = 0;
                        else
                            U[i, j] = (Matrix[i, j] - sum) / U[i, i];
                    }
                }
            }
            var UT = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    UT[i, j] = U[j, i];
                }
            }


            var flag = true;
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    if (double.IsNaN(U[j, k]) || double.IsNaN(UT[j, k])) flag = false;
                }
            }
            if (flag)
            {
                var x = SolveSlau(UT, U, B);
                Console.WriteLine("Решение системы:");
                PrintResult(x);
            }
            else
            {
                Console.WriteLine("Не существует UTU разложения");
            }
        }



        //СТЕПЕННОЙ МЕТОД
        public void PowerIteration()
        {
            var invMatrix = GetInvMatrix();
            var eigenvector = new double[N];
            var eigenValue = GetEigenvalue(invMatrix, eigenvector);
            Console.WriteLine("Минимальное собственное число матрицы: {0:F3}", eigenValue);
            Console.WriteLine("Cоотвествующий собственный вектор:");
            PrintResult(eigenvector);
        }

        //получение минимального собственного числа и соотвествующего вектора
        private double GetEigenvalue(double[,] matrix, double[] eigenvector)
        {
            Console.WriteLine("Задайте точность вычислений для степенного метода:");
            var eps = double.Parse(Console.ReadLine());
            var X = new double[N];
            var X0 = new double[N];
            var X0norm = new double[N];
            var d = 0.0;
            var d0 = 0.0;
            var e = 1.0;
            X0[0] = 1;
            for (int i = 1; i < N; i++)
            {
                X0[i] = 0;
            }
            do
            {
                var sum = 0.0;
                for (int i = 0; i < N; i++)
                {
                    sum += X0[i] * X0[i];
                }
                d0 = Math.Sqrt(sum);
                for (int i = 0; i < N; i++)
                {
                    X0norm[i] = X0[i] / d0;
                }
                for (int i = 0; i < N; i++)
                {
                    X[i] = 0;
                    for (int j = 0; j < N; j++)
                    {
                        X[i] += matrix[i, j] * X0norm[j];
                    }
                }
                sum = 0;
                for (int i = 0; i < N; i++)
                {
                    sum += X[i] * X[i];
                }
                d = Math.Sqrt(sum);
                e = Math.Abs(d - d0);
                for (int i = 0; i < N; i++)
                {
                    X0[i] = X[i];
                }
            } while (e > eps);
            X0norm = Rationing(X0norm);
            for (int i = 0; i < N; i++)
            {
                eigenvector[i] = X0norm[i];
            }
            return 1/d;
        }

        //получение обратной матрицы
        private double[,] GetInvMatrix()
        {
            var matrix = Matrix;
            var determinant = GetDeterminant(matrix);
            var invMatrix = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    var minor = GetMinor(i, j, matrix);
                    invMatrix[j, i] = Math.Pow(-1, i + j)*GetDeterminant(minor) / determinant;
                }
            }
            return invMatrix;
        }

        //получение миноров
        private double[,] GetMinor(int i, int j, double[,] matrix)
        {
            var minor = new double[N - 1, N - 1];
            var rowRemoved = 0;
            var colRemoved = 0;
            for (int k = 0; k < matrix.GetLength(0) - 1; k++)
            {
                if (k == i)
                {
                    rowRemoved = 1;
                }
                colRemoved = 0;
                for (int l = 0; l < matrix.GetLength(0) - 1; l++)
                {

                    if (l == j)
                    {
                        colRemoved = 1;
                    }
                    minor[k, l] = matrix[k + rowRemoved, l + colRemoved];
                }
            }
            return minor;
        }

        //определитель матрицы
        private double GetDeterminant(double[,] matrix)
        {
            var L = new double[matrix.GetLength(0), matrix.GetLength(0)];
            var U = new double[matrix.GetLength(0), matrix.GetLength(0)];
            GetLU(L, U, matrix.GetLength(0), matrix);
            var determinant = 1.0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                determinant *= L[i, i];
                determinant *= U[i, i];
            }
            if (determinant == 0) throw new Exception("Матрица вырожденная");
            return determinant;
        }

        //нормирование
        private double[] Rationing(double[] eigenvector)
        {
            var sum = 0.0;
            for (int i = 0; i < N; i++)
            {
                sum += eigenvector[i] * eigenvector[i];
            }
            for (int i = 0; i < N; i++)
            {
                eigenvector[i] = eigenvector[i] / Math.Sqrt(sum);
            }
            return eigenvector;
        }



        //ПОЛЯ КЛАССА
        //коэффициенты 
        public double[,] Matrix { get; set; }

        //массив свободных членов
        public double[] B { get; set; }

        //размерность матрицы
        public int N { get; set; }
    }
}
