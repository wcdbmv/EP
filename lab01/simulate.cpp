#include "simulate.hpp"

#include <algorithm>
#include <functional>
#include <list>
#include <memory>
#include <numeric>
#include <queue>
#include <vector>

#include "random.hpp"


class IConsumer {
public:
	virtual ~IConsumer() = default;

	virtual bool Receive(double current_event_time) = 0;
};

class IProducer {
public:
	virtual ~IProducer() = default;

	virtual void Attach(std::shared_ptr<IConsumer> consumer) = 0;
	virtual bool Process(double current_event_time) = 0;
};


class DummyConsumer : public IConsumer {
public:
	bool Receive(double /* current_event_time */) override {
		throw std::runtime_error("DummyConsumer->Receive(double) should never be called");
		return true;
	}
};

class DummyProducer : public IProducer {
public:
	void Attach(std::shared_ptr<IConsumer> consumer) override {
		consumers_.push_back(consumer);
	}

	bool Process(double /* current_event_time */) override {
		return true;
	}

protected:
	std::list<std::shared_ptr<IConsumer>> consumers_;
};


class RequestGenerator : public DummyProducer, public DummyConsumer {
public:
	using Distribution = std::function<double()>;

	explicit RequestGenerator(Distribution&& distribution)
		: distribution_(std::move(distribution))
		, next_event_time_(0.0)
		, n_generated_(0)
	{
	}

	[[maybe_unused]] virtual double GenerateNextEventTime(double current_event_time) {
		++n_generated_;
		return next_event_time_ = current_event_time + distribution_();
	}

	bool Process(double current_event_time) override {
		GenerateNextEventTime(current_event_time);
		for (auto&& consumer : consumers_) {

			consumer->Receive(current_event_time);
			return true;
		}
		return false;
	}

	[[nodiscard]] double GetNextEventTime() const {
		return next_event_time_;
	}

	[[nodiscard]] size_t GetNGenerated() const {
		return n_generated_;
	}

protected:
	Distribution distribution_;
	double next_event_time_;
	size_t n_generated_;
};

class RequestProcessor : public RequestGenerator {
public:
	explicit RequestProcessor(Distribution&& distribution, size_t queue_size_limit = 0)
		: RequestGenerator(std::forward<Distribution>(distribution))
		, queue_size_limit_(queue_size_limit)
		, n_processed_(0)
		, queue_()
		, waiting_times_()
		, loading_times_()
	{
	}

	bool Process(double current_event_time) override {
		const auto cond = !queue_.empty();
		if (cond) {
			waiting_times_.push_back(current_event_time - queue_.front());
			queue_.pop();

			++n_processed_;
			RequestGenerator::Process(current_event_time);
		}
		return cond;
	}

	bool Receive(double current_event_time) override {
		const auto cond = queue_size_limit_ == 0 || queue_.size() < queue_size_limit_;
		if (cond) {
			queue_.push(current_event_time);
			RequestGenerator::GenerateNextEventTime(current_event_time);
			loading_times_.push_back(next_event_time_ - current_event_time);
		}
		return cond;
	}

	[[maybe_unused]] double GenerateNextEventTime(double current_event_time) override {
		next_event_time_ = queue_.empty() ? 0.0 : RequestGenerator::GenerateNextEventTime(current_event_time);
		if (next_event_time_ > 0.0) {
			loading_times_.push_back(next_event_time_ - current_event_time);
		}
		return next_event_time_;
	}

	[[nodiscard]] size_t GetNProcessed() const {
		return n_processed_;
	}

	[[nodiscard]] const std::vector<double>& GetWaitingTimes() const {
		return waiting_times_;
	}

	[[nodiscard]] const std::vector<double>& GetLoadingTimes() const {
		return loading_times_;
	}

protected:
	size_t queue_size_limit_;
	size_t n_processed_;
	std::queue<double> queue_;
	std::vector<double> waiting_times_;
	std::vector<double> loading_times_;
};


class Simulator {
public:
	[[maybe_unused]] std::shared_ptr<RequestGenerator> AddDevice(std::shared_ptr<RequestGenerator> device) {
		devices_.push_back(device);
		return device;
	}

	void Attach(size_t producer, size_t consumer) {
		if (producer == consumer || producer >= devices_.size() || consumer >= devices_.size()) {
			throw std::invalid_argument("Simulator.Attach(): invalid indices");
		}

		devices_[producer]->Attach(devices_[consumer]);
	}

	std::shared_ptr<RequestGenerator> GetDevice(size_t index) {
		return devices_[index];
	}

	void Simulate(double until) {
		for (auto current_event_time = GetCurrentEventTime(); current_event_time <= until; current_event_time = GetCurrentEventTime()) {
			for (auto&& device : devices_) {
				if (std::abs(current_event_time - device->GetNextEventTime()) < 1e-12) {
					device->Process(current_event_time);
				}
			}
		}
	}

private:
	std::vector<std::shared_ptr<RequestGenerator>> devices_;

	double GetCurrentEventTime() const {
		return std::accumulate(devices_.begin(), devices_.end(), 9999999999.0, [](auto&& acc, auto& run) {
			return 0.0 < run->GetNextEventTime() && run->GetNextEventTime() < acc ? run->GetNextEventTime() : acc;
		});
	}
};


SimulateResult Simulate(const SimulateParams& params) {
	Simulator simulator;
	simulator.AddDevice(std::make_shared<RequestGenerator>([=]() { return uniform_real(params.a, params.b); }))->GenerateNextEventTime(0.0);
	simulator.AddDevice(std::make_shared<RequestProcessor>([=]() { return weibull_real(params.k, params.lambda); }));
	simulator.Attach(0, 1);
	simulator.Simulate(params.t);

	auto processor = std::dynamic_pointer_cast<RequestProcessor>(simulator.GetDevice(1));
	const auto& loading_times = processor->GetLoadingTimes();
	const auto loading_time_of_processor = std::accumulate(loading_times.begin(), loading_times.end(), 0.0);
	return {
		.load = std::min(loading_time_of_processor / params.t, 1.0),
	};
}
