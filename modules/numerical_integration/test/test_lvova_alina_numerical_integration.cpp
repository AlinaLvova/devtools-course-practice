// Copyright 2018 Lvova Alina

#include <gtest/gtest.h>
#include <string>
#include "include/numerical_integration.h"

double norm(double x) {
    return (x * x * x / 3 + x * x / 2 + 2 * x);
}

TEST(Lvova_Alina_NumericalIntegrationTest, CorrectDivisionsInConstructor) {
    // Arrange
    int  N = 100;

    // Act & Assert
    ASSERT_NO_THROW(Integral I(3.04, 9.45, N));
}

TEST(Lvova_Alina_NumericalIntegrationTest, IncorrectLimitsInConstructor) {
    // Arrange
    double A = 3.14, B = -65.654;

    // Act & Assert
    try {
        Integral I(A, B, 10);
    }
    catch (std::out_of_range const & err) {
        EXPECT_EQ(err.what(), std::string("limits out of range"));
    }
}

TEST(Lvova_Alina_NumericalIntegrationTest, IncorrectLimitInConstructor) {
    // Arrange
    double A = 3.14, B = 65.654;

    // Act & Assert
    ASSERT_NO_THROW(Integral I(3.04, 9.45, 100));
}

TEST(Lvova_Alina_NumericalIntegrationTest, RiemannSumLeft) {
    // Arrange
    double A = -13.14, B = 65.654;
    int  N = 1000;
    double result, norma;
    Integral I(A, B, N);
    double epsilon = 1.8e2;

    // Act
    result = I.RiemannSumLeft();
    norma = norm(B) - norm(A);

    // Assert
    EXPECT_NEAR(norma, result, epsilon);
}

TEST(Lvova_Alina_NumericalIntegrationTest, TrapezoidalRule) {
    // Arrange
    double A = -13.14, B = 65.654;
    int  N = 1000;
    double result, norma;
    Integral I(A, B, N);
    double epsilon = 82e-3;

    // Act
    result = I.TrapezoidalRule();
    norma = norm(B) - norm(A);

    // Assert
    EXPECT_NEAR(norma, result, epsilon);
}

TEST(Lvova_Alina_NumericalIntegrationTest, SimpsonRule) {
    // Arrange
    double A = -153.26485, B = 35.365695;
    int N = 750;
    double result, norma;
    Integral I(A, B, N);
    double epsilon = 3.26e-9;

    // Act
    result = I.SimpsonRule();
    norma = norm(B) - norm(A);

    // Assert
    EXPECT_NEAR(norma, result, epsilon);
}

TEST(Lvova_Alina_NumericalIntegrationTest, Simpson3_8Rule) {
    // Arrange
    double A = -13.14, B = 65.654;
    int  N = 1000;
    double result, norma;
    Integral I(A, B, N);
    double epsilon = 2.65e-9;

    // Act
    result = I.Simpson3_8Rule();
    norma = norm(B) - norm(A);

    // Assert
    EXPECT_NEAR(norma, result, epsilon);
}

TEST(Lvova_Alina_NumericalIntegrationTest, BooleRule) {
    // Arrange
    double A = -13.14, B = 65.654;
    int  N = 1000;
    double result, norma;
    Integral I(A, B, N);
    double epsilon = 2.67e-9;

    // Act
    result = I.BooleRule();
    norma = norm(B) - norm(A);

    // Assert
    EXPECT_NEAR(norma, result, epsilon);
}

TEST(Lvova_Alina_NumericalIntegrationTest, NewtonCotes5) {
    // Arrange
    double A = -13.14, B = 65.654;
    int  N = 1000;
    double result, norma;
    Integral I(A, B, N);
    double epsilon = 2.65e-9;

    // Act
    result = I.NewtonCotes5();
    norma = norm(B) - norm(A);

    // Assert
    EXPECT_NEAR(norma, result, epsilon);
}

TEST(Lvova_Alina_NumericalIntegrationTest, SetCorrectLowerAndCorrectDivisions) {
    // Arrange
    double A = -153.26485, B = 0;
    int N = 750;
    double result, norma;
    Integral I;
    double epsilon = 3.5e-9;

    // Act
    I.setLower(A);
    I.setDivisions(N);
    result = I.GaussianQuadrature();
    norma = norm(B) - norm(A);

    // Assert
    EXPECT_NEAR(norma, result, epsilon);
}

TEST(Lvova_Alina_NumericalIntegrationTest, SetIncorrectLower) {
    // Arrange
    double A = -153.26485, B = 10;
    int N = 750;
    Integral I(A, B, N);
    double newLower = 10.10;

    // Act & Assert
    try {
        I.setLower(newLower);
    }
    catch (std::out_of_range const & err) {
        EXPECT_EQ(err.what(), std::string("lower limit out of range"));
    }
}

TEST(Lvova_Alina_NumericalIntegrationTest, SetCorrectUpperAndCorrectDivisions) {
    // Arrange
    double A = 0, B = 35.365695;
    int N = 750;
    double result, norma;
    Integral I;
    double epsilon = 3.5e-9;

    // Act
    I.setUpper(B);
    I.setDivisions(N);
    result = I.GaussianQuadrature();
    norma = norm(B) - norm(A);

    // Assert
    EXPECT_NEAR(norma, result, epsilon);
}

TEST(Lvova_Alina_NumericalIntegrationTest, SetIncorrectUpper) {
    // Arrange
    double A = -15.25, B = 10;
    int N = 750;
    Integral I(A, B, N);
    double newUpper = -28.10;

    // Act & Assert
    try {
        I.setUpper(newUpper);
    }
    catch (std::out_of_range const & err) {
        EXPECT_EQ(err.what(), std::string("upper limit out of range"));
    }
}

TEST(Lvova_Alina_NumericalIntegrationTest, SetIncorrectDivisions) {
    // Arrange
    Integral I;
    double newDiv = -1;

    // Act & Assert
    try {
        I.setDivisions(newDiv);
    }
    catch (std::out_of_range const & err) {
        EXPECT_EQ(err.what(), std::string("Divisions must be greater than 1"));
    }
}
