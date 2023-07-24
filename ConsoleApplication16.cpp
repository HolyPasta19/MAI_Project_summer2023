#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <cctype> // Для использования std::isdigit

class ComplexNumber {
public:
    double real;
    double imaginary;

    ComplexNumber(double real = 0.0, double imaginary = 0.0) : real(real), imaginary(imaginary) {}

    ComplexNumber operator+(const ComplexNumber& other) const {
        return ComplexNumber(real + other.real, imaginary + other.imaginary);
    }

    ComplexNumber operator-(const ComplexNumber& other) const {
        return ComplexNumber(real - other.real, imaginary - other.imaginary);
    }

    ComplexNumber operator*(const ComplexNumber& other) const {
        double realPart = real * other.real - imaginary * other.imaginary;
        double imagPart = real * other.imaginary + imaginary * other.real;
        return ComplexNumber(realPart, imagPart);
    }
};
bool isMatrixSquare(const std::vector<std::vector<ComplexNumber>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    return rows == cols;
}
std::vector<std::vector<ComplexNumber>> enterMatrixManually(int& rows, int& cols) {
    std::cout << "Введите количество строк матрицы: ";
    std::cin >> rows;
    std::cout << "Введите количество столбцов матрицы: ";
    std::cin >> cols;

    if (rows <= 0 || cols <= 0) {
        std::cerr << "Некорректный размер матрицы. Размеры должны быть положительными числами." << std::endl;
        return enterMatrixManually(rows, cols); // Ввод матрицы вручную
    }

    if (rows != cols) {
        std::cerr << "Некорректный размер матрицы. Матрица должна быть квадратной (количество строк равно количеству столбцов)." << std::endl;
        return enterMatrixManually(rows, cols); // Ввод матрицы вручную
    }

    std::vector<std::vector<ComplexNumber>> matrix(rows, std::vector<ComplexNumber>(cols));

    std::cout << "Введите элементы матрицы:\n";
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double realPart, imagPart;
            std::cout << "Элемент (" << i + 1 << ", " << j + 1 << "): ";
            std::cin >> realPart >> imagPart;
            matrix[i][j] = ComplexNumber(realPart, imagPart);
        }
    }

    return matrix;
}

std::vector<std::vector<ComplexNumber>> readMatrixFromFile(const std::string& filename, int& rows, int& cols) {
    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Ошибка при открытии файла " << filename << std::endl;
        return enterMatrixManually(rows, cols); // Ввод матрицы вручную
    }

    inputFile >> rows >> cols;

    if (rows <= 0 || cols <= 0) {
        std::cerr << "Некорректные размеры матрицы в файле " << filename << std::endl;
        inputFile.close();
        return enterMatrixManually(rows, cols); // Ввод матрицы вручную
    }

    if (rows != cols) {
        std::cerr << "Некорректный размер матрицы в файле " << filename << ". Матрица должна быть квадратной (количество строк равно количеству столбцов)." << std::endl;
        inputFile.close();
        return enterMatrixManually(rows, cols); // Ввод матрицы вручную
    }

    std::vector<std::vector<ComplexNumber>> matrix(rows, std::vector<ComplexNumber>(cols));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double realPart, imagPart;
            char plusSign, iSymbol;

            if (!(inputFile >> realPart)) {
                std::cerr << "Ошибка в чтении данных из файла: " << filename << std::endl;
                inputFile.close();
                return enterMatrixManually(rows, cols); // Ввод матрицы вручную
            }

            inputFile >> plusSign;
            if (plusSign != '+') {
                std::cerr << "Ошибка в формате комплексного числа: " << filename << std::endl;
                inputFile.close();
                return enterMatrixManually(rows, cols); // Ввод матрицы вручную
            }

            inputFile >> imagPart;
            inputFile >> iSymbol;

            if (iSymbol != 'i') {
                std::cerr << "Ошибка в формате комплексного числа: " << filename << std::endl;
                inputFile.close();
                return enterMatrixManually(rows, cols); // Ввод матрицы вручную
            }

            matrix[i][j] = ComplexNumber(realPart, imagPart);
        }
    }

    inputFile.close();
    return matrix;
}



int getNearestPowerOfTwo(int num) {
    int power = 1;
    while (power < num) {
        power *= 2;
    }
    return power;
}
bool isNotPowerOfTwo(int n) {
    return !(n > 0 && (n & (n - 1)) == 0);
}

std::vector<std::vector<ComplexNumber>> addMatrix(const std::vector<std::vector<ComplexNumber>>& A, const std::vector<std::vector<ComplexNumber>>& B) {
    int n = A.size();
    std::vector<std::vector<ComplexNumber>> C(n, std::vector<ComplexNumber>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<ComplexNumber>> subtractMatrix(const std::vector<std::vector<ComplexNumber>>& A, const std::vector<std::vector<ComplexNumber>>& B) {
    int n = A.size();
    std::vector<std::vector<ComplexNumber>> C(n, std::vector<ComplexNumber>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<ComplexNumber>> strassenMultiplication(const std::vector<std::vector<ComplexNumber>>& A, const std::vector<std::vector<ComplexNumber>>& B) {
    int n = A.size();

    if (n <= 2) {
        std::vector<std::vector<ComplexNumber>> C(n, std::vector<ComplexNumber>(n, 0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    C[i][j] = C[i][j] + A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }

    int newSize = n / 2;

    std::vector<std::vector<ComplexNumber>> A11(newSize, std::vector<ComplexNumber>(newSize));
    std::vector<std::vector<ComplexNumber>> A12(newSize, std::vector<ComplexNumber>(newSize));
    std::vector<std::vector<ComplexNumber>> A21(newSize, std::vector<ComplexNumber>(newSize));
    std::vector<std::vector<ComplexNumber>> A22(newSize, std::vector<ComplexNumber>(newSize));

    std::vector<std::vector<ComplexNumber>> B11(newSize, std::vector<ComplexNumber>(newSize));
    std::vector<std::vector<ComplexNumber>> B12(newSize, std::vector<ComplexNumber>(newSize));
    std::vector<std::vector<ComplexNumber>> B21(newSize, std::vector<ComplexNumber>(newSize));
    std::vector<std::vector<ComplexNumber>> B22(newSize, std::vector<ComplexNumber>(newSize));

    for (int i = 0; i < newSize; ++i) {
        for (int j = 0; j < newSize; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + newSize];
            A21[i][j] = A[i + newSize][j];
            A22[i][j] = A[i + newSize][j + newSize];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + newSize];
            B21[i][j] = B[i + newSize][j];
            B22[i][j] = B[i + newSize][j + newSize];
        }
    }

    std::vector<std::vector<ComplexNumber>> M1 = strassenMultiplication(addMatrix(A11, A22), addMatrix(B11, B22));
    std::vector<std::vector<ComplexNumber>> M2 = strassenMultiplication(addMatrix(A21, A22), B11);
    std::vector<std::vector<ComplexNumber>> M3 = strassenMultiplication(A11, subtractMatrix(B12, B22));
    std::vector<std::vector<ComplexNumber>> M4 = strassenMultiplication(A22, subtractMatrix(B21, B11));
    std::vector<std::vector<ComplexNumber>> M5 = strassenMultiplication(addMatrix(A11, A12), B22);
    std::vector<std::vector<ComplexNumber>> M6 = strassenMultiplication(subtractMatrix(A21, A11), addMatrix(B11, B12));
    std::vector<std::vector<ComplexNumber>> M7 = strassenMultiplication(subtractMatrix(A12, A22), addMatrix(B21, B22));

    std::vector<std::vector<ComplexNumber>> C11 = subtractMatrix(addMatrix(addMatrix(M1, M4), M7), M5);
    std::vector<std::vector<ComplexNumber>> C12 = addMatrix(M3, M5);
    std::vector<std::vector<ComplexNumber>> C21 = addMatrix(M2, M4);
    std::vector<std::vector<ComplexNumber>> C22 = subtractMatrix(subtractMatrix(addMatrix(M1, M3), M2), M6);

    std::vector<std::vector<ComplexNumber>> C(n, std::vector<ComplexNumber>(n, 0));
    for (int i = 0; i < newSize; ++i) {
        for (int j = 0; j < newSize; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + newSize] = C12[i][j];
            C[i + newSize][j] = C21[i][j];
            C[i + newSize][j + newSize] = C22[i][j];
        }
    }

    return C;
}

std::vector<std::vector<ComplexNumber>> padMatrixWithZeros(const std::vector<std::vector<ComplexNumber>>& matrix, int newSize) {

    int oldRows = matrix.size();
    int oldCols = matrix[0].size();

    std::vector<std::vector<ComplexNumber>> newMatrix(newSize, std::vector<ComplexNumber>(newSize, 0));

    for (int i = 0; i < oldRows; ++i) {
        for (int j = 0; j < oldCols; ++j) {
            newMatrix[i][j] = matrix[i][j];
        }
    }

    return newMatrix;
}

#include <complex>


void printMatrix(const std::vector<std::vector<ComplexNumber>>& matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << matrix[i][j].real << " + " << matrix[i][j].imaginary << "i\t";
        }
        std::cout << std::endl;
    }
}
int main() {
    setlocale(LC_ALL, "Russian");

    int rowsA, colsA, rowsB, colsB;

    std::cout << "Матрица A:\n";
    std::vector<std::vector<ComplexNumber>> matrixA;

    int choiceA;
    std::cout << "Выберите способ ввода матрицы A:\n";
    std::cout << "1. Ввести матрицу A вручную.\n";
    std::cout << "2. Прочитать матрицу A из файла.\n";
    std::cout << "Ваш выбор (1 или 2): ";
    std::cin >> choiceA;

    switch (choiceA) {
    case 1:
        matrixA = enterMatrixManually(rowsA, colsA);
        break;
    case 2: {
        std::string filenameA;
        std::cout << "Введите имя файла для матрицы A: ";
        std::cin >> filenameA;
        matrixA = readMatrixFromFile(filenameA, rowsA, colsA);
        break;
    }
    default:
        std::cout << "Некорректный выбор. Выбран ввод матрицы A вручную по умолчанию.\n";
        matrixA = enterMatrixManually(rowsA, colsA);
        break;
    }
    std::cout << "Матрица A:\n";
    printMatrix(matrixA);


    std::cout << "Матрица B:\n";
    std::vector<std::vector<ComplexNumber>> matrixB;

    int choiceB;
    std::cout << "Выберите способ ввода матрицы B:\n";
    std::cout << "1. Ввести матрицу B вручную.\n";
    std::cout << "2. Прочитать матрицу B из файла.\n";
    std::cout << "Ваш выбор (1 или 2): ";
    std::cin >> choiceB;

    switch (choiceB) {
    case 1:
        matrixB = enterMatrixManually(rowsB, colsB);
        break;
    case 2: {
        std::string filenameB;
        std::cout << "Введите имя файла для матрицы B: ";
        std::cin >> filenameB;
        matrixB = readMatrixFromFile(filenameB, rowsB, colsB);
        break;
    }
    default:
        std::cout << "Некорректный выбор. Выбран ввод матрицы B вручную по умолчанию.\n";
        matrixB = enterMatrixManually(rowsB, colsB);
        break;
    }
    std::cout << "Матрица B:\n";
    printMatrix(matrixB);

    // Проверяем, возможно ли умножить матрицы A и B
    if (colsA != rowsB) {
        std::cout << "Умножение матриц A и B невозможно из-за несовпадения размерностей." << std::endl;
        return 1;
    }
    if (rowsA != colsA || rowsB != colsB || isNotPowerOfTwo(rowsA) || isNotPowerOfTwo(rowsB)) {
        int choice;
        std::cout << "Матрицы имеют некорректный формат, что может сказаться на скорости работы программы.\n";
        std::cout << "Матрицы будут дополнены нулями.\n";
        std::cout << "Желаете продолжить?:\n";
        std::cout << "1. Да\n";
        std::cout << "2. Нет\n";
        std::cout << "Ваш выбор (1 или 2): ";
        std::cin >> choice;
        if (choice != 1) {
            std::cout << "Программа завершила свою работу." << std::endl;
            return 0;
        }

        int newSize = std::pow(2, std::max(std::ceil(std::log2(rowsA)), std::ceil(std::log2(colsB))));
        matrixA = padMatrixWithZeros(matrixA, newSize);
        matrixB = padMatrixWithZeros(matrixB, newSize);
    }

    // Умножаем матрицы методом Штрассена
    std::vector<std::vector<ComplexNumber>> resultMatrix = strassenMultiplication(matrixA, matrixB);

    // Выводим результат умножения
    std::cout << "Результат умножения матриц A и B:\n";
    printMatrix(resultMatrix);

    // Вводим имя файла для сохранения результата
    std::string resultFilename;
    std::cout << "Введите имя файла для сохранения результата: ";
    std::cin >> resultFilename;

    // Сохраняем результат в файл
    std::ofstream outputFile(resultFilename);
    if (!outputFile) {
        std::cerr << "Ошибка при создании файла " << resultFilename << std::endl;
        return 1;
    }

    outputFile << resultMatrix.size() << " " << resultMatrix[0].size() << std::endl;
    for (const auto& row : resultMatrix) {
        for (const auto& num : row) {
            outputFile << num.real << " + " << num.imaginary << "i ";
        }
        outputFile << std::endl;
    }

    outputFile.close();

    return 0;
}
