#include "simulate.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <list>
#include <memory>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <vector>

#include "random.hpp"
#include "statistics.hpp"


struct Event {
	double time;
	uint8_t type;
};

bool operator==(const Event& lhs, const Event& rhs) {
	return std::equal_to<double>{}(lhs.time, rhs.time) && lhs.type == rhs.type;
}


class IConsumer {
public:
	virtual ~IConsumer() = default;

	virtual bool Receive(const Event& next_event) = 0;
};

class IProducer {
public:
	virtual ~IProducer() = default;

	virtual void Attach(std::shared_ptr<IConsumer> consumer) = 0;
	virtual bool Process(const Event& current_event) = 0;
};


class DummyConsumer : public IConsumer {
public:
	bool Receive(const Event& /* next_event */) override {
		throw std::runtime_error("DummyConsumer->Receive(const Event&) should never be called");
		return true;
	}
};

class DummyProducer : public IProducer {
public:
	void Attach(std::shared_ptr<IConsumer> consumer) override {
		consumers_.push_back(consumer);
	}

	bool Process(const Event& /* current_event */) override {
		throw std::runtime_error("DummyProducer->Process(const Event&) should never be called");
		return true;
	}

protected:
	std::list<std::shared_ptr<IConsumer>> consumers_;
};


class IDevice : public DummyProducer, public DummyConsumer {
public:
	using Distribution = std::function<double()>;

	[[nodiscard]] virtual Event GetNextEvent() const = 0;
};


class RequestGenerator : public IDevice {
public:
	using Distribution = IDevice::Distribution;

	explicit RequestGenerator(Distribution&& distribution, uint8_t type)
		: distribution_{std::move(distribution)}
		, next_event_time_{distribution_()}
		, type_{type}
	{
	}

	// current_event.time == next_event_time_
	bool Process(const Event& current_event) override {
		assert(current_event == GetNextEvent());

		GenerateNextEventTime(current_event.time);
		for (auto&& consumer : consumers_) {
			if (consumer->Receive(current_event)) {
				return true;
			}
		}
		return false;
	}

	[[nodiscard]] Event GetNextEvent() const override {
		return {next_event_time_, type_};
	}

private:
	Distribution distribution_;
	double next_event_time_;
	uint8_t type_;

	void GenerateNextEventTime(double current_event_time) {
		next_event_time_ = current_event_time + distribution_();
	}
};


class RequestProcessor : public IDevice {
public:
	using Distribution = IDevice::Distribution;
	using Distributions = std::unordered_map<uint8_t, Distribution>;

	RequestProcessor(Distributions&& distributions)
		: distributions_{std::move(distributions)}
		, next_event_{}
		, queue_{}
		, waiting_times_{}
		, n_processed_{}
	{
	}

	bool Receive(const Event& current_event) override {
		queue_.push(current_event);
		if (std::equal_to<double>{}(next_event_.time, 0.0)) {
			GenerateNextEventTimeOnReceive(current_event);
		}
		return true;
	}

	// current_event.time == next_event_time_
	bool Process(const Event& current_event) override {
		assert(current_event == next_event_);
		assert(!queue_.empty());

		waiting_times_.push_back(current_event.time - queue_.front().time);
		queue_.pop();

		++n_processed_;
		GenerateNextEventTimeOnProcess(current_event);

		return true;
	}

	[[nodiscard]] Event GetNextEvent() const override {
		return next_event_;
	}

	[[nodiscard]] size_t GetNProcessed() const {
		return n_processed_;
	}

	[[nodiscard]] double GetAverageWaitingTime() const {
		return std::accumulate(waiting_times_.begin(), waiting_times_.end(), 0.0) / static_cast<double>(waiting_times_.size());
	}

private:
	Distributions distributions_;
	Event next_event_;
	std::queue<Event> queue_;
	std::vector<double> waiting_times_;
	size_t n_processed_;

	void GenerateNextEventTimeOnReceive(const Event& current_event) {
		next_event_ = {current_event.time + distributions_[current_event.type](), current_event.type};
	}

	void GenerateNextEventTimeOnProcess(const Event& current_event) {
		const auto may_process = !queue_.empty();
		if (may_process) {
			GenerateNextEventTimeOnReceive(current_event);
		} else {
			next_event_ = Event{};
		}
	}
};


class Simulator {
public:
	[[maybe_unused]] std::shared_ptr<IDevice> AddDevice(std::shared_ptr<IDevice> device) {
		devices_.push_back(device);
		return device;
	}

	void Attach(size_t producer, size_t consumer) {
		if (producer == consumer || producer >= devices_.size() || consumer >= devices_.size()) {
			throw std::invalid_argument("Simulator.Attach(): invalid indices");
		}

		devices_[producer]->Attach(devices_[consumer]);
	}

	std::shared_ptr<IDevice> GetDevice(size_t index) {
		return devices_[index];
	}

	void Simulate(double until) {
		for (auto current_event = FindCurrentEvent(); current_event.time <= until; current_event = FindCurrentEvent()) {
			for (auto&& device : devices_) {
				if (current_event == device->GetNextEvent()) {
					device->Process(current_event);
					break;
				}
			}
		}
	}

private:
	std::vector<std::shared_ptr<IDevice>> devices_;

	Event FindCurrentEvent() const {
		return std::accumulate(devices_.begin(), devices_.end(), Event{9999999999.0, 0}, [](auto&& acc, auto& run) {
			const auto event = run->GetNextEvent();
			const auto time = event.time;
			return 0.0 < time && time < acc.time ? event : acc;
		});
	}
};


SimulateResult Simulate(const SimulateParams& params) {
	Simulator simulator;
	simulator.AddDevice(std::make_shared<RequestGenerator>([=]() { return uniform_real(params.a1, params.b1); }, 1));
	simulator.AddDevice(std::make_shared<RequestGenerator>([=]() { return uniform_real(params.a2, params.b2); }, 2));
	simulator.AddDevice(std::make_shared<RequestProcessor>(RequestProcessor::Distributions{
		{1, [=]() { return weibull_real(params.k1, params.lambda1); }},
		{2, [=]() { return weibull_real(params.k2, params.lambda2); }},
	}));
	simulator.Attach(0, 2);
	simulator.Attach(1, 2);
	simulator.Simulate(params.t);

	auto processor = std::dynamic_pointer_cast<RequestProcessor>(simulator.GetDevice(2));
	return {
		.average_waiting = processor->GetAverageWaitingTime(),
	};
}
