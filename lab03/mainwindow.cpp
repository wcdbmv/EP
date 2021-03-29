#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include "random.hpp"


MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	const auto setUpGroupBox = [](const auto& range, auto* minLineEdit, auto* maxLineEdit, auto* doubleSpinBox) {
		minLineEdit->setText(QString::number(range.min));
		maxLineEdit->setText(QString::number(range.max));
		// doubleSpinBox->setMinimum(min);
		// doubleSpinBox->setMaximum(max);
		doubleSpinBox->setValue((range.max + range.min) / 2.0);
	};

	setUpGroupBox(FFE_PARAMS.lambda1, ui->lambda1MinLineEdit, ui->lambda1MaxLineEdit, ui->lambda1DoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.lambda2, ui->lambda2MinLineEdit, ui->lambda2MaxLineEdit, ui->lambda2DoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.mu, ui->muMinLineEdit, ui->muMaxLineEdit, ui->muDoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.sigma_lambda1, ui->sigmaLambda1MinLineEdit, ui->sigmaLambda1MaxLineEdit, ui->sigmaLambda1DoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.sigma_lambda2, ui->sigmaLambda2MinLineEdit, ui->sigmaLambda2MaxLineEdit, ui->sigmaLambda2DoubleSpinBox);

	const auto setUpTableWidget = [](auto* tableWidget) {
		tableWidget->verticalHeader()->setVisible(false);
		tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
	};

	setUpTableWidget(ui->fullFactorialExperimentTableWidget);
	setUpTableWidget(ui->fractionalFactorialExperimentTableWidget);

	ui->fullFactorialExperimentTableWidget->verticalHeader()->setVisible(false);
	ui->fullFactorialExperimentTableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

	full_factorial_experiment_result = FullFactorialExperiment(FFE_PARAMS);
	fractional_factorial_experiment_result = FractionalFactorialExperiment(FFE_PARAMS);

	const auto insertRow = [](auto* tableWidget, const FfeTableRow& row) {
		const auto rows = tableWidget->rowCount();
		tableWidget->insertRow(rows);
		tableWidget->setItem(rows, 0, new QTableWidgetItem(QString::number(row.index)));
		tableWidget->setItem(rows, 1, new QTableWidgetItem(QString::number(row.x1)));
		tableWidget->setItem(rows, 2, new QTableWidgetItem(QString::number(row.x2)));
		tableWidget->setItem(rows, 3, new QTableWidgetItem(QString::number(row.x3)));
		tableWidget->setItem(rows, 4, new QTableWidgetItem(QString::number(row.x4)));
		tableWidget->setItem(rows, 5, new QTableWidgetItem(QString::number(row.x5)));
		tableWidget->setItem(rows, 6, new QTableWidgetItem(QString::number(row.y_mean)));
		tableWidget->setItem(rows, 7, new QTableWidgetItem(QString::number(row.y_var)));
		tableWidget->setItem(rows, 8, new QTableWidgetItem(QString::number(row.y_hat)));
		tableWidget->setItem(rows, 9, new QTableWidgetItem(QString::number(row.dy_hat)));
		tableWidget->setItem(rows, 10, new QTableWidgetItem(QString::number(row.u_hat)));
		tableWidget->setItem(rows, 11, new QTableWidgetItem(QString::number(row.du_hat)));
	};

	for (auto&& row : full_factorial_experiment_result.table) {
		insertRow(ui->fullFactorialExperimentTableWidget, row);
	}

	for (auto&& row : fractional_factorial_experiment_result.table) {
		insertRow(ui->fractionalFactorialExperimentTableWidget, row);
	}

	const auto setUpRegressionLineEdit = [](auto* lineEdit, const PartialNonlinearCoefficients<5>& cf, size_t limit = 1) {
		constexpr size_t K = 5;
		assert(limit <= K);

		const auto coef = [](double coefficient) {
			if (coefficient < 0) {
				return " - " + QString::number(-coefficient);
			}
			return " + " + QString::number(coefficient);
		};

		const auto x = [](size_t i) {
			return QString{"*x"} + QString::number(i);
		};

		QString s = QString::number(cf[0]);
		for (size_t i = 0; (limit >= 1) && i < 5; ++i) {

		}

		// hell no
		for (size_t k = 1; k < K; ++k) {
			for (size_t a = 1; limit >= 1 && a <= K; ++a) {
				if (k == 1) {
					s += coef(cf.a(a)) + x(a);
				} else {
					for (size_t b = a + 1; limit >= 2 && b <= K; ++b) {
						if (k == 2) {
							s += coef(cf.a(a, b)) + x(a) + x(b);
						} else {
							for (size_t c = b + 1; limit >= 3 && c <= K; ++c) {
								if (k == 3) {
									s += coef(cf.a(a, b, c)) + x(a) + x(b) + x(c);
								} else {
									for (size_t d = c + 1; limit >= 4 && d <= K; ++d) {
										s += coef(cf.a(a, b, c, d)) + x(a) + x(b) + x(c) + x(d);
									}
								}
							}
						}
					}
				}

			}
		}
		if (limit >= 5) {
			s += coef(cf.a(1, 2, 3, 4, 5)) + x(1) + x(2) + x(3) + x(4) + x(5);
		}

		lineEdit->setText(s);
	};

	setUpRegressionLineEdit(ui->yHatRegressionFullFactorialExperimentLineEdit, full_factorial_experiment_result.coefficients, 1);
	setUpRegressionLineEdit(ui->uHatRegressionFullFactorialExperimentLineEdit, full_factorial_experiment_result.coefficients, 5);
	setUpRegressionLineEdit(ui->yHatRegressionFractionalFactorialExperimentLineEdit, fractional_factorial_experiment_result.coefficients, 1);
	setUpRegressionLineEdit(ui->uHatRegressionFractionalFactorialExperimentLineEdit, fractional_factorial_experiment_result.coefficients, 5);
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::on_calculatePushButton_clicked()
{
	const DotParameters dot_params{
		.lambda1       = ui->lambda1DoubleSpinBox->value(),
		.lambda2       = ui->lambda2DoubleSpinBox->value(),
		.mu            = ui->muDoubleSpinBox->value(),
		.sigma_lambda1 = ui->sigmaLambda1DoubleSpinBox->value(),
		.sigma_lambda2 = ui->sigmaLambda2DoubleSpinBox->value(),
	};

	const auto actual = CalculateDot(dot_params, FFE_PARAMS.times);

	const auto factors = NormalizeFactors(FFE_PARAMS, dot_params);
	const auto y_hat_full = CalculateDotWithRegression(full_factorial_experiment_result.coefficients, factors, 1);
	const auto u_hat_full = CalculateDotWithRegression(full_factorial_experiment_result.coefficients, factors, 5);
	const auto y_hat_frac = CalculateDotWithRegression(fractional_factorial_experiment_result.coefficients, factors, 1);
	const auto u_hat_frac = CalculateDotWithRegression(fractional_factorial_experiment_result.coefficients, factors, 5);

	ui->actualAverageWaitingTimeLineEdit->setText(QString::number(actual));
	ui->yHatFullFactorialExperimentLineEdit->setText(QString::number(y_hat_full));
	ui->uHatFullFactorialExperimentLineEdit->setText(QString::number(u_hat_full));
	ui->yHatFractionalFactorialExperimentLineEdit->setText(QString::number(y_hat_frac));
	ui->uHatFractionalFactorialExperimentLineEdit->setText(QString::number(u_hat_frac));
	ui->dyHatFullFactorialExperimentLineEdit->setText(QString::number(actual - y_hat_full));
	ui->duHatFullFactorialExperimentLineEdit->setText(QString::number(actual - u_hat_full));
	ui->dyHatFractionalFactorialExperimentLineEdit->setText(QString::number(actual - y_hat_frac));
	ui->duHatFractionalFactorialExperimentLineEdit->setText(QString::number(actual - u_hat_frac));
}
