#!/usr/bin/env bash
set -euo pipefail

# PETSC_OPTIONS may also contain solver overrides, for example:
# PETSC_OPTIONS="-ksp_type cg -pc_type hypre -pc_hypre_type boomeramg"
# Add -mpm_petsc_warm_start to benchmark nonlinear-iteration warm starts.
# Add -ksp_monitor there to print per-iteration residuals.
# MPM_PROFILE_MAX_TIME and MPM_PROFILE_MAX_TIMESTEPS can override the short
# default simulation limits; MPM_PROFILE_TIMEOUT_SECONDS overrides the wall
# timeout. The input file is copied to a temporary file.
if (($# == 0)); then
  echo "Usage: $0 <command> [arguments...]" >&2
  echo "Example: $0 ../../dbg_petsc_improve/bin/vaango input.ups" >&2
  exit 2
fi

profile_max_time="${MPM_PROFILE_MAX_TIME:-2.e-8}"
profile_max_timesteps="${MPM_PROFILE_MAX_TIMESTEPS:-2}"
profile_timeout="${MPM_PROFILE_TIMEOUT_SECONDS:-60}"
if [[ ! "$profile_max_time" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
  echo "MPM_PROFILE_MAX_TIME must be a positive numeric value." >&2
  exit 2
fi
if [[ ! "$profile_max_timesteps" =~ ^[1-9][0-9]*$ ]]; then
  echo "MPM_PROFILE_MAX_TIMESTEPS must be a positive integer." >&2
  exit 2
fi
if [[ ! "$profile_timeout" =~ ^[1-9][0-9]*$ ]]; then
  echo "MPM_PROFILE_TIMEOUT_SECONDS must be a positive integer." >&2
  exit 2
fi
if ! command -v timeout >/dev/null 2>&1; then
  echo "The 'timeout' command is required to bound the profile run." >&2
  exit 127
fi

profile_options="-mpm_petsc_profile -ksp_converged_reason"
if [[ -n "${PETSC_OPTIONS:-}" ]]; then
  export PETSC_OPTIONS="${PETSC_OPTIONS} ${profile_options}"
else
  export PETSC_OPTIONS="${profile_options}"
fi

profile_args=("$@")
profile_input=""
for ((i = 1; i < ${#profile_args[@]}; i++)); do
  input_arg="${profile_args[$i]}"
  if [[ "$input_arg" == *.ups && -f "$input_arg" ]]; then
    input_dir="$(cd -- "$(dirname -- "$input_arg")" && pwd -P)"
    input_file="${input_dir}/$(basename -- "$input_arg")"
    if [[ $(grep -c '<maxTime>' "$input_file" || true) -ne 1 ]]; then
      echo "Expected exactly one <maxTime> entry in $input_file." >&2
      exit 2
    fi

    profile_input="$(mktemp "${input_dir}/.profile-implicit-XXXXXX.ups")"
    awk -v max_time="$profile_max_time" -v max_steps="$profile_max_timesteps" '
      /<maxTime>/ {
        sub(/<maxTime>[^<]*<\/maxTime>/, "<maxTime>" max_time "</maxTime>")
      }
      /<max_timesteps>/ {
        sub(/<max_timesteps>[^<]*<\/max_timesteps>/,
            "<max_timesteps>" max_steps "</max_timesteps>")
        has_max_steps = 1
      }
      /<\/Time>/ {
        if (!has_max_steps) {
          print "       <max_timesteps>" max_steps "</max_timesteps>"
          has_max_steps = 1
        }
      }
      { print }
    ' "$input_file" >"$profile_input"

    profile_args[$i]="$profile_input"
    trap 'rm -f -- "$profile_input"' EXIT
    break
  fi
done

if [[ -z "$profile_input" ]]; then
  echo "Warning: no .ups input found; simulation time limits were not overridden." >&2
fi

timeout --signal=TERM --kill-after=5s "${profile_timeout}s" "${profile_args[@]}"
