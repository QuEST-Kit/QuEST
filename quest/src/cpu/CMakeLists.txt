# @author Oliver Thomson Brown

target_sources(QuEST
  PRIVATE
  cpu_config.cpp
  cpu_subroutines.cpp
)

# issue #598: optionally use compensated (Kahan) summation in applyCompMatr's
# dense-matrix CPU subroutine. Off by default (naive summation); enable with
#   -DQUEST_COMPENSATE_DENSEMATR_SUM=ON
# to trade ~2-7x runtime for ~25x better accuracy on many-target matrices.
option(QUEST_COMPENSATE_DENSEMATR_SUM "Use compensated (Kahan) summation in applyCompMatr's CPU subroutine (issue #598)" OFF)
if (QUEST_COMPENSATE_DENSEMATR_SUM)
  target_compile_definitions(QuEST PRIVATE QUEST_COMPENSATED_DENSEMATR_SUM=1)
  message(STATUS "Compensated (Kahan) dense-matrix summation is turned ON.")
endif()