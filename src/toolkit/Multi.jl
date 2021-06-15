# Multi
Spar <: Array{ }
println(typeof(Multi))
println(dump(Multi))

Multi.parameters = Core.svec(5,5,5,5,5)

#=
DataType <: Type{T}
  name::Core.TypeName
  super::DataType
  parameters::Core.SimpleVector
  types::Core.SimpleVector
  names::Core.SimpleVector
  instance::Any
  layout::Ptr{Nothing}
  size::Int32
  ninitialized::Int32
  hash::Int32
  abstract::Bool
  mutable::Bool
  hasfreetypevars::Bool
  isconcretetype::Bool
  isdispatchtuple::Bool
  isbitstype::Bool
  zeroinit::Bool
  isinlinealloc::Bool
  has_concrete_subtype::Bool
=#


type SparseArray{Tv,Ti,N} <: AbstractSparseArray{Tv,Ti,N} where {Ti <: NTuple{N,<:Integer}, Tv, N}
