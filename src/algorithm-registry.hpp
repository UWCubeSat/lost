#ifndef ALGORITHM_REGISTRY_H
#define ALGORITHM_REGISTRY_H

#include <string>
#include <map>
#include <functional>
#include <memory>

namespace lost {

class PipelineOptions;

/**
 * A registry that maps string names to factory functions for a given algorithm base class.
 * Each algorithm category (centroiding, star-id, attitude estimation) gets its own registry
 * instance. Algorithms register themselves at static-initialization time in their .cpp files,
 * so adding a new algorithm never requires modifying SetPipeline().
 */
template <typename T>
class AlgorithmRegistry {
public:
    using FactoryFn = std::function<std::unique_ptr<T>(const PipelineOptions &)>;

    static AlgorithmRegistry &Instance() {
        static AlgorithmRegistry instance;
        return instance;
    }

    void Register(const std::string &name, FactoryFn factory) {
        factories_[name] = std::move(factory);
    }

    /// Create an algorithm by name, or return nullptr if the name is not registered.
    std::unique_ptr<T> Create(const std::string &name, const PipelineOptions &options) const {
        auto it = factories_.find(name);
        if (it == factories_.end()) {
            return nullptr;
        }
        return it->second(options);
    }

    bool Has(const std::string &name) const {
        return factories_.count(name) > 0;
    }

private:
    AlgorithmRegistry() = default;
    std::map<std::string, FactoryFn> factories_;
};

} // namespace lost

#endif
