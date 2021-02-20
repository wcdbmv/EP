#include "simulate.hpp"

#include <algorithm>
#include <functional>
#include <list>
#include <memory>
#include <numeric>
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
		return next_event_time_ = current_event_time + distribution_();
	}

	bool Process(double current_event_time) override {
		++n_generated_;
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
	explicit RequestProcessor(Distribution&& distribution, size_t max_queue_size = 0.0)
		: RequestGenerator(std::forward<Distribution>(distribution))
		, max_queue_size_(max_queue_size)
		, cur_queue_size_(0)
		, n_processed_(0)
	{
	}

	bool Process(double current_event_time) override {
		const auto cond = cur_queue_size_ > 0;
		if (cond) {
			--cur_queue_size_;
			++n_processed_;
			RequestGenerator::Process(current_event_time);
		}
		return cond;
	}

	bool Receive(double current_event_time) override {
		const auto cond = max_queue_size_ == 0 || cur_queue_size_ < max_queue_size_;
		if (cond) {
			++cur_queue_size_;
			RequestGenerator::GenerateNextEventTime(current_event_time);
		}
		return cond;
	}

	[[maybe_unused]] double GenerateNextEventTime(double current_event_time) override {
		return next_event_time_
			= cur_queue_size_
			? RequestGenerator::GenerateNextEventTime(current_event_time)
			: 0.0;
	}

protected:
	size_t max_queue_size_;
	size_t cur_queue_size_;
	size_t n_processed_;
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
		double current_event_time = 0.0;
		while (current_event_time < until) {
			current_event_time = std::accumulate(devices_.begin(), devices_.end(), 0.0, [](auto&& acc, auto& run) {
				return 0.0 < run->GetNextEventTime() && run->GetNextEventTime() < acc ? run->GetNextEventTime() : acc;
			});

			for (auto&& device : devices_) {
				if (std::abs(current_event_time - device->GetNextEventTime()) < 1e-12) {
					device->Process(current_event_time);
				}
			}
		}
	}

private:
	std::vector<std::shared_ptr<RequestGenerator>> devices_;
};


SimulateResult Simulate(const SimulateParams& params) {
	Simulator simulator;
	simulator.AddDevice(std::make_shared<RequestGenerator>([=]() { return uniform_real(params.a, params.b); }))->GenerateNextEventTime(0.0);
	simulator.AddDevice(std::make_shared<RequestGenerator>([=]() { return weibull_real(params.k, params.lambda); }));
	simulator.Attach(0, 1);
	simulator.Simulate(params.t);
	return {};
}
