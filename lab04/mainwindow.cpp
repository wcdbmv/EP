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

	setUpGroupBox(EXP_PARAMS.lambda1, ui->lambda1MinLineEdit, ui->lambda1MaxLineEdit, ui->lambda1DoubleSpinBox);
	setUpGroupBox(EXP_PARAMS.lambda2, ui->lambda2MinLineEdit, ui->lambda2MaxLineEdit, ui->lambda2DoubleSpinBox);
	setUpGroupBox(EXP_PARAMS.mu1, ui->mu1MinLineEdit, ui->mu1MaxLineEdit, ui->mu1DoubleSpinBox);
	setUpGroupBox(EXP_PARAMS.mu2, ui->mu2MinLineEdit, ui->mu2MaxLineEdit, ui->mu2DoubleSpinBox);
	setUpGroupBox(EXP_PARAMS.sigma_lambda1, ui->sigmaLambda1MinLineEdit, ui->sigmaLambda1MaxLineEdit, ui->sigmaLambda1DoubleSpinBox);
	setUpGroupBox(EXP_PARAMS.sigma_lambda2, ui->sigmaLambda2MinLineEdit, ui->sigmaLambda2MaxLineEdit, ui->sigmaLambda2DoubleSpinBox);

	exp_result = RunExperiment(EXP_PARAMS);

	const auto setUpTableWidget = [&](auto* tableWidget) {
		tableWidget->verticalHeader()->setVisible(false);
		tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

		tableWidget->setColumnCount(static_cast<int>(exp_result.planning_matrix[0].size()));

		constexpr size_t K = 6;
		const auto x = [](size_t i) { return QString{"x"} + QString::number(i); };
		const auto x2mS = [&](size_t i) { return x(i) + "^2 - S"; };
		QStringList labels{"№"};
		for (size_t i = 0; i <= 6; ++i) {
			labels << x(i);
		}
		for (auto&& combination : GenerateAllCombinations(K, 2)) {
			labels << x(combination[0] + 1) + x(combination[1] + 1);
		}
		for (size_t i = 1; i <= 6; ++i) {
			labels << x2mS(i);
		}
		labels << "y" << "ŷ" << "|y - ŷ|";

		tableWidget->setHorizontalHeaderLabels(labels);
	};

	setUpTableWidget(ui->occdTableWidget);

	const auto insertRow = [](auto* tableWidget, const std::vector<double>& row) {
		const auto rows = tableWidget->rowCount();
		tableWidget->insertRow(rows);
		for (size_t i = 0; i < row.size(); ++i) {
			tableWidget->setItem(rows, static_cast<int>(i), new QTableWidgetItem(QString::number(row[i])));
		}
	};

	for (auto&& row : exp_result.planning_matrix) {
		insertRow(ui->occdTableWidget, row);
	}

	const auto setUpRegressionLineEdit = [](auto* lineEdit, const std::vector<double>& cf) {
		constexpr size_t K = 6;

		const auto coef = [](double coefficient) {
			return coefficient < 0 ? " - " + QString::number(-coefficient) : " + " + QString::number(coefficient);
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

		const double S = std::sqrt(64.0 / 77.0);

		const auto xmS = [S](size_t index) {
			return QString{"*(x"} + QString::number(index + 1) + QString{"² - "} + QString::number(S) + QString{")"};
		};

		for (size_t i = 0; i < K; ++i) {
			double c = qAbs(cf[j++]);
			c = i / 2 == 1 ? -c : c;
			s += coef(c) + xmS(i);
		}

		lineEdit->setText(s);
	};

	setUpRegressionLineEdit(ui->yHatRegressionOccdLineEdit, exp_result.coefficients);
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

	const auto actual = CalculateDot(dot_params, EXP_PARAMS.times);

	const auto factors = NormalizeFactors(EXP_PARAMS, dot_params);
	const auto y_hat_full = CalculateDotWithRegression(exp_result.coefficients, factors);

	ui->actualAverageWaitingTimeLineEdit->setText(QString::number(actual));
	ui->yHatOccdLineEdit->setText(QString::number(y_hat_full));
	ui->dyHatOccdLineEdit->setText(QString::number(qAbs(actual - y_hat_full)));
}
