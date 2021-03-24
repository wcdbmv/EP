#include "simulate.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <list>
#include <memory>
#include <numeric>
#include <queue>
#include <vector>

#include "random.hpp"
#include "statistics.hpp"


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
		throw std::runtime_error("DummyProducer->Process(double) should never be called");
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

	// current_event_time == next_event_time_
	bool Process(double current_event_time) override {
		GenerateNextEventTime(current_event_time);
		for (auto&& consumer : consumers_) {
			if (consumer->Receive(current_event_time)) {
				return true;
			}
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

	// current_event_time == next_event_time_
	bool Process(double current_event_time) override {
		const auto can_pop = !queue_.empty();
		if (can_pop) {
			waiting_times_.push_back(current_event_time - queue_.front());
			queue_.pop();

			++n_processed_;
			RequestGenerator::Process(current_event_time);
		}
		return can_pop;
	}

	bool Receive(double current_event_time) override {
		const auto can_push = !HasLimit() || queue_.size() < queue_size_limit_;
		if (can_push) {
			queue_.push(current_event_time);
			if (next_event_time_ <= 0.0) {
				RequestGenerator::GenerateNextEventTime(current_event_time);
				loading_times_.push_back(next_event_time_ - current_event_time);
			}
		}
		return can_push;
	}

	[[maybe_unused]] double GenerateNextEventTime(double current_event_time) override {
		const auto may_process = !queue_.empty();
		if (may_process) {
			RequestGenerator::GenerateNextEventTime(current_event_time);
			loading_times_.push_back(next_event_time_ - current_event_time);
		} else {
			next_event_time_ = 0.0;
		}
		return next_event_time_;
	}

	[[nodiscard]] size_t GetNProcessed() const {
		return n_processed_;
	}

	[[nodiscard]] size_t GetQueueSize() const {
		return queue_.size();
	}

	[[nodiscard]] const std::vector<double>& GetWaitingTimes() const {
		return waiting_times_;
	}

	[[nodiscard]] double GetAverageWaitingTime() const {
		return std::accumulate(waiting_times_.begin(), waiting_times_.end(), 0.0) / static_cast<double>(waiting_times_.size());
	}

	[[nodiscard]] const std::vector<double>& GetLoadingTimes() const {
		return loading_times_;
	}

	[[nodiscard]] double GetLoadTime() const {
		return std::accumulate(loading_times_.begin(), loading_times_.end(), 0.0);
	}

	[[nodiscard]] double GetAverageLoadTime() const {
		return GetLoadTime() / static_cast<double>(loading_times_.size());
	}

	// this is estimate actual load time (xD)
	[[nodiscard]] double GetExtrapolatedLoadTime() const {
		return GetLoadTime() + GetAverageLoadTime() * static_cast<double>(queue_.size());
	}

protected:
	size_t queue_size_limit_;
	size_t n_processed_;
	std::queue<double> queue_;
	std::vector<double> waiting_times_;
	std::vector<double> loading_times_;

	[[nodiscard]] bool HasLimit() const {
		return queue_size_limit_;
	}
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
		for (auto current_event_time = FindCurrentEventTime(); current_event_time <= until; current_event_time = FindCurrentEventTime()) {
			for (auto&& device : devices_) {
				if (std::abs(current_event_time - device->GetNextEventTime()) < 1e-12) {
					device->Process(current_event_time);
				}
			}
		}
	}

private:
	std::vector<std::shared_ptr<RequestGenerator>> devices_;

	double FindCurrentEventTime() const {
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
	return {
		.average_waiting = processor->GetAverageWaitingTime(),
	};
}
