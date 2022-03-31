const DEFAULT_ENCODINGS = [
    :domain_wall,
    :one_hot,
]

const DEFAULT_OPTIMIZERS = [
    :descent,
]

const DEFAULT_LEARNING_RATES = [
    0.1,
    0.01,
    0.001,
]

const DEFAULT_DOMAIN_SIZES = [
    3,
    5,
    10,
]

const DEFAULT_TRAINING_SET_SIZES = [
    100,
    1000,
    10000,
    100000,
]

const DEFAULT_ROUNDIND_PRECISIONS = [
    100,
    10,
    5,
    1,
]

const ALL_PARAMETERS = Dict(
    # Search space parameters
    :domain_size => DEFAULT_DOMAIN_SIZES,
    :encoding => DEFAULT_ENCODINGS,
    :training_size_set => DEFAULT_TRAINING_SET_SIZES,

    # Learning parameters
    :learning_rate => DEFAULT_LEARNING_RATES,
    :optimizer => DEFAULT_OPTIMIZERS,
    :rouding_precision => DEFAULT_ROUNDIND_PRECISIONS,
)
