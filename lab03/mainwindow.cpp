#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include "random.hpp"
#include <cstdlib>


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
	setUpGroupBox(FFE_PARAMS.mu1, ui->mu1MinLineEdit, ui->mu1MaxLineEdit, ui->mu1DoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.mu2, ui->mu2MinLineEdit, ui->mu2MaxLineEdit, ui->mu2DoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.sigma_lambda1, ui->sigmaLambda1MinLineEdit, ui->sigmaLambda1MaxLineEdit, ui->sigmaLambda1DoubleSpinBox);
	setUpGroupBox(FFE_PARAMS.sigma_lambda2, ui->sigmaLambda2MinLineEdit, ui->sigmaLambda2MaxLineEdit, ui->sigmaLambda2DoubleSpinBox);

	const auto setUpTableWidget = [](auto* tableWidget) {
		tableWidget->verticalHeader()->setVisible(false);
		tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
	};

	setUpTableWidget(ui->fullFactorialExperimentTableWidget);
	setUpTableWidget(ui->fractionalFactorialExperimentTableWidget);

	full_result = FullFactorialExperiment(FFE_PARAMS);
	frac_result = FractionalFactorialExperiment(FFE_PARAMS);

	const auto insertRow = [](auto* tableWidget, const FfeTableRow& row) {
		const auto rows = tableWidget->rowCount();
		tableWidget->insertRow(rows);
		tableWidget->setItem(rows, 0, new QTableWidgetItem(QString::number(row.index)));
		tableWidget->setItem(rows, 1, new QTableWidgetItem(QString::number(row.x1)));
		tableWidget->setItem(rows, 2, new QTableWidgetItem(QString::number(row.x2)));
		tableWidget->setItem(rows, 3, new QTableWidgetItem(QString::number(row.x3)));
		tableWidget->setItem(rows, 4, new QTableWidgetItem(QString::number(row.x4)));
		tableWidget->setItem(rows, 5, new QTableWidgetItem(QString::number(row.x5)));
		tableWidget->setItem(rows, 6, new QTableWidgetItem(QString::number(row.x6)));
		tableWidget->setItem(rows, 7, new QTableWidgetItem(QString::number(row.y_mean)));
		tableWidget->setItem(rows, 8, new QTableWidgetItem(QString::number(row.y_var)));
		tableWidget->setItem(rows, 9, new QTableWidgetItem(QString::number(row.y_hat)));
		tableWidget->setItem(rows, 10, new QTableWidgetItem(QString::number(row.dy_hat)));
		tableWidget->setItem(rows, 11, new QTableWidgetItem(QString::number(row.u_hat)));
		tableWidget->setItem(rows, 12, new QTableWidgetItem(QString::number(row.du_hat)));
	};

	for (auto&& row : full_result.table) {
		insertRow(ui->fullFactorialExperimentTableWidget, row);
	}

//	const auto insertRow2 = [](auto* tableWidget, const FfeTableRow& row) {
//		const auto rows = tableWidget->rowCount();
//		volatile int sign = (-1) + 2 * (rand() < 0.5);
//		tableWidget->insertRow(rows);
//		tableWidget->setItem(rows, 0, new QTableWidgetItem(QString::number(row.index)));
//		tableWidget->setItem(rows, 1, new QTableWidgetItem(QString::number(row.x1)));
//		tableWidget->setItem(rows, 2, new QTableWidgetItem(QString::number(row.x2)));
//		tableWidget->setItem(rows, 3, new QTableWidgetItem(QString::number(row.x3)));
//		tableWidget->setItem(rows, 4, new QTableWidgetItem(QString::number(row.x4)));
//		tableWidget->setItem(rows, 5, new QTableWidgetItem(QString::number(row.x5)));
//		tableWidget->setItem(rows, 6, new QTableWidgetItem(QString::number(row.x6)));
//		tableWidget->setItem(rows, 7, new QTableWidgetItem(QString::number(row.y_mean)));
//		tableWidget->setItem(rows, 8, new QTableWidgetItem(QString::number(row.y_var)));
//		tableWidget->setItem(rows, 9, new QTableWidgetItem(QString::number(row.y_hat)));
//		tableWidget->setItem(rows, 10, new QTableWidgetItem(QString::number(row.dy_hat)));
//		tableWidget->setItem(rows, 11, new QTableWidgetItem(QString::number(row.y_mean + row.y_mean * sign * 0.03)));
//		tableWidget->setItem(rows, 12, new QTableWidgetItem(QString::number(row.y_mean * sign * 0.03)));
//	};

	for (auto&& row : frac_result.table) {
		insertRow(ui->fractionalFactorialExperimentTableWidget, row);
	}

	const auto setUpRegressionLineEdit = [](auto* lineEdit, const PartialNonlinearCoefficients<6>& cf, size_t limit) {
		constexpr size_t K = 6;
		assert(limit <= K);

		const auto coef = [](double coefficient) {
			if (coefficient < 0) {
				return " - " + QString::number(-coefficient);
			}
			return " + " + QString::number(coefficient);
		};

		const auto xs = [](const std::vector<size_t>& combination) {
			QString result = "";
			for (auto index : combination) {
				result += QString{"*x"} + QString::number(index + 1);
			}
			return result;
		};

		QString s = QString::number(cf[0]);

		for (size_t i = 1, j = 1; i <= limit; ++i) {
			for (auto&& combination : GenerateAllCombinations(K, i)) {
				s += coef(cf[j++]) + xs(combination);
			}
		}

		lineEdit->setText(s);
	};

	setUpRegressionLineEdit(ui->yHatRegressionFullFactorialExperimentLineEdit, full_result.coefficients, 1);
	setUpRegressionLineEdit(ui->uHatRegressionFullFactorialExperimentLineEdit, full_result.coefficients, 6);
	setUpRegressionLineEdit(ui->yHatRegressionFractionalFactorialExperimentLineEdit, frac_result.coefficients, 1);
	setUpRegressionLineEdit(ui->uHatRegressionFractionalFactorialExperimentLineEdit, frac_result.coefficients, 2);
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
		.mu1           = ui->mu1DoubleSpinBox->value(),
		.mu2           = ui->mu2DoubleSpinBox->value(),
		.sigma_lambda1 = ui->sigmaLambda1DoubleSpinBox->value(),
		.sigma_lambda2 = ui->sigmaLambda2DoubleSpinBox->value(),
	};

	const auto actual = CalculateDot(dot_params, FFE_PARAMS.times);

	const auto factors = NormalizeFactors(FFE_PARAMS, dot_params);
	const auto y_hat_full = CalculateDotWithRegression(full_result.coefficients, factors, 1);
	const auto u_hat_full = CalculateDotWithRegression(full_result.coefficients, factors, 6);
	const auto y_hat_frac = CalculateDotWithRegression(frac_result.coefficients, factors, 1);
	const auto u_hat_frac = CalculateDotWithRegression(frac_result.coefficients, factors, 2);

	ui->actualAverageWaitingTimeLineEdit->setText(QString::number(actual));
	ui->yHatFullFactorialExperimentLineEdit->setText(QString::number(y_hat_full));
	ui->uHatFullFactorialExperimentLineEdit->setText(QString::number(u_hat_full));
	ui->yHatFractionalFactorialExperimentLineEdit->setText(QString::number(y_hat_frac));
	ui->uHatFractionalFactorialExperimentLineEdit->setText(QString::number(u_hat_frac));
	ui->dyHatFullFactorialExperimentLineEdit->setText(QString::number(qAbs(actual - y_hat_full)));
	ui->duHatFullFactorialExperimentLineEdit->setText(QString::number(qAbs(actual - u_hat_full)));
	ui->dyHatFractionalFactorialExperimentLineEdit->setText(QString::number(qAbs(actual - y_hat_frac)));
	ui->duHatFractionalFactorialExperimentLineEdit->setText(QString::number(qAbs(actual - u_hat_frac)));
}
