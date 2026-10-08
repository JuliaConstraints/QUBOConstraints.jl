module QUBOConstraints

# SECTION - Usings and imports
using ConstraintDomains
using LinearAlgebra
using TestItems
using SparseArrays
import SHA
import TOML

# SECTION - exports
export binarize
export debinarize
export is_valid
export train

export QUBO_linear_sum
export BitID, Codebook, DecodeResult, codebook, encode_code, encoded_codes, decode_code
export QUBOComponent, energy, dense_qubo, canonical_equal, compose, rename_auxiliaries
export ExhaustiveReport, exhaustive_check, constraint_oracle, component_artifact, write_component
export one_hot_constraint
export one_hot_validity, guard_one_hot
export AbstractCodebook, StructuredCodebook, structured_codebook, encoding_registry, encoding_validity
export UnaryRelation, BinaryRelation, local_component, local_constraint, table_constraint
export AnyValue, ANY_VALUE, guard_encodings
export truth_component, CompilationLimit
export IntegerCondition, satisfies, sum_component, count_component
export nvalues_component, cardinality_component, extremum_component, knapsack_component
export binpacking_component, finite_arithmetic_component
export table_component, regular_component, mdd_component, slide_scopes, slide_component
export lex_component, element_component, precedence_component, stretch_component
export cumulative_component, nooverlap_component, circuit_component
export IntensionNode, intension_component
export AtomicNode, AtomicPlan, AtomicCompilation, atomic_plan, compile_atomic, atomic_values, atomic_witness
export AtomicCondition, core_atomic_expression, core_atomic_plan
export validate_atomic_plan
export AtomicSearchSpace, atomic_decision_domains, atomic_operator_code
export encode_atomic_weights, decode_atomic_weights, compile_atomic_weights
export AtomicEnumerativeOptimizer, AtomicTrainingResult
export ValuePairGuidance

# SECTION - includes
include("base.jl")
include("codebook.jl")
include("structured_codebook.jl")
include("component.jl")
include("oracle.jl")
include("artifacts.jl")
include("one_hot_constraints.jl")
include("validity.jl")
include("encoding_validity.jl")
include("truth_component.jl")
include("local_relations.jl")
include("arithmetic_components.jl")
include("aggregate_components.jl")
include("structural_components.jl")
include("finite_structural.jl")
include("scheduling_components.jl")
include("intension_components.jl")
include("atomic_definedness.jl")
include("atomic_components.jl")
include("atomic_recipes.jl")

include("handmade/linear_sum.jl")

include("encoding/domain_wall.jl")
include("encoding/one_hot.jl")
include("encoding/conversion.jl")

include("learn.jl")
include("atomic_search.jl")
include("value_pair_guidance.jl")

end
