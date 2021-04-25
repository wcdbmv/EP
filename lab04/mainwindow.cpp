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

	setUpGroupBox(OCCD_PARAMS.lambda1, ui->lambda1MinLineEdit, ui->lambda1MaxLineEdit, ui->lambda1DoubleSpinBox);
	setUpGroupBox(OCCD_PARAMS.lambda2, ui->lambda2MinLineEdit, ui->lambda2MaxLineEdit, ui->lambda2DoubleSpinBox);
	setUpGroupBox(OCCD_PARAMS.mu1, ui->mu1MinLineEdit, ui->mu1MaxLineEdit, ui->mu1DoubleSpinBox);
	setUpGroupBox(OCCD_PARAMS.mu2, ui->mu2MinLineEdit, ui->mu2MaxLineEdit, ui->mu2DoubleSpinBox);
	setUpGroupBox(OCCD_PARAMS.sigma_lambda1, ui->sigmaLambda1MinLineEdit, ui->sigmaLambda1MaxLineEdit, ui->sigmaLambda1DoubleSpinBox);
	setUpGroupBox(OCCD_PARAMS.sigma_lambda2, ui->sigmaLambda2MinLineEdit, ui->sigmaLambda2MaxLineEdit, ui->sigmaLambda2DoubleSpinBox);

	const auto setUpTableWidget = [](auto* tableWidget) {
		tableWidget->verticalHeader()->setVisible(false);
		tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
	};

	setUpTableWidget(ui->occdTableWidget);

	occd_result = OrthogonalCentralCompositeDesign(OCCD_PARAMS);;

	const auto insertRow = [](auto* tableWidget, const OccdTableRow& row) {
		const auto rows = tableWidget->rowCount();
		tableWidget->insertRow(rows);
		tableWidget->setItem(rows, 0, new QTableWidgetItem(QString::number(row.index)));
		tableWidget->setItem(rows, 1, new QTableWidgetItem(QString::number(row.x1)));
		tableWidget->setItem(rows, 2, new QTableWidgetItem(QString::number(row.x2)));
		tableWidget->setItem(rows, 3, new QTableWidgetItem(QString::number(row.x3)));
		tableWidget->setItem(rows, 4, new QTableWidgetItem(QString::number(row.x4)));
		tableWidget->setItem(rows, 5, new QTableWidgetItem(QString::number(row.x5)));
		tableWidget->setItem(rows, 6, new QTableWidgetItem(QString::number(row.x6)));
		tableWidget->setItem(rows, 7, new QTableWidgetItem(QString::number(row.x12mS)));
		tableWidget->setItem(rows, 8, new QTableWidgetItem(QString::number(row.x22mS)));
		tableWidget->setItem(rows, 9, new QTableWidgetItem(QString::number(row.x32mS)));
		tableWidget->setItem(rows, 10, new QTableWidgetItem(QString::number(row.x42mS)));
		tableWidget->setItem(rows, 11, new QTableWidgetItem(QString::number(row.x52mS)));
		tableWidget->setItem(rows, 12, new QTableWidgetItem(QString::number(row.x62mS)));
		tableWidget->setItem(rows, 13, new QTableWidgetItem(QString::number(row.y_mean)));
		tableWidget->setItem(rows, 14, new QTableWidgetItem(QString::number(row.y_var)));
		tableWidget->setItem(rows, 15, new QTableWidgetItem(QString::number(row.y_hat)));
		tableWidget->setItem(rows, 16, new QTableWidgetItem(QString::number(row.dy_hat)));
	};

	for (auto&& row : occd_result.table) {
		insertRow(ui->occdTableWidget, row);
	}

	const auto setUpRegressionLineEdit = [](auto* lineEdit, const NonlinearCoefficients<6>& cf) {
		constexpr size_t K = 6;

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

		size_t j = 1;
		for (size_t i = 1; i <= 2; ++i) {
			for (auto&& combination : GenerateAllCombinations(K, i)) {
				s += coef(cf[j++]) + xs(combination);
			}
		}

		const double S = std::sqrt(6.0 / 77.0);

		const auto xmS = [S](size_t index) {
			return QString{"*(x"} + QString::number(index + 1) + QString{"^2 - "} + QString::number(S) + QString{")"};
		};

		for (size_t i = 0; i < K; ++i) {
			s += coef(cf[j++]) + xmS(i);
		}

		lineEdit->setText(s);
	};

	setUpRegressionLineEdit(ui->yHatRegressionOccdLineEdit, occd_result.coefficients);
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

	const auto actual = CalculateDot(dot_params, OCCD_PARAMS.times);

	const auto factors = NormalizeFactors(OCCD_PARAMS, dot_params);
	const auto y_hat_full = CalculateDotWithRegression(occd_result.coefficients, factors);

	ui->actualAverageWaitingTimeLineEdit->setText(QString::number(actual));
	ui->yHatOccdLineEdit->setText(QString::number(y_hat_full));
	ui->dyHatOccdLineEdit->setText(QString::number(qAbs(actual - y_hat_full)));
}
