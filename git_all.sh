#!/bin/bash

# Array of submodule paths
paths=(
"CMakeConfigs"
"SBB/AutoCorr"
"SBB/AutoCorr/aCorrsOTF"
"SBB/AutoCorr/aCorrsOTF/CMakeConfigs"
"SBB/FFT/FFTW_extra"
"SBB/FFT/FFTW_extra/CMakeConfigs"
"SBB/FFT/FFT_CUDA"
"SBB/FFT/FFT_CUDA/CMakeConfigs"
"SBB/Histograms"
"SBB/Histograms/CMakeConfigs"
"SBB/Math_extra"
"SBB/Math_extra/CMakeConfigs"
"SBB/Multi_array"
"SBB/Multi_array/CMakeConfigs"
"SBB/Omp_extra"
"SBB/Omp_extra/CMakeConfigs"
"SBB/Time_quadratures"
"SBB/Time_quadratures/CMakeConfigs"
)

#Sort directories deepest first.
IFS=$'\n' paths=($(sort -r <<<"${paths[*]}"))
unset IFS
#printf "%s\n" "${sorted[@]}"

#Docstring
usage() {
 echo "Usage: $0 [OPTIONS]"
 echo "Options:"
 echo " -h, --help  , help      Display this help message"
 echo " -r, --rebase, rebase    Rebase all submodule to a given branch"
 echo " -c, --commit, commit    Commits all submodule and the parent module using the given message"
}

# Check if not options given
if [ $# -eq 0 ]; then
    echo "No arguments provided. Please provide either 'rebase' or 'commit'."
	echo ""
	usage
    exit 1
fi

function special_print() {
    local line="#####"
    local length=${#1}
    for ((i=1; i<$length; i++)); do
        line+="#"
    done
    echo "$line"
    echo "# $1 #"
    echo "$line"
}

function checkout_branch() {
    local branch="$1"
    if git rev-parse --verify "$branch" >/dev/null 2>&1; then
        git checkout "$branch"
    else
        echo "Error: Branch '$branch' does not exist."
        exit 1
    fi
}

function checkout_paths() {
    local branch="$1"
    local paths=("${@:2}")
	special_print "Checkout to $branch"
    for path in "${paths[@]}"; do
		echo "$path"
        git -C "$path" checkout_branch "$branch"
    done
}

function pull_paths() {
    local paths=("$@")
    special_print "Pull"
    for path in "${paths[@]}"; do
        echo "$path"
        git -C "$path" pull
    done
}

function push_paths() {
    local paths=("$@")
    special_print "Push"
    for path in "${paths[@]}"; do
        echo "$path"
        git -C "$path" push
    done
}

function force_push_paths() {
    local paths=("$@")
    special_print "Force Push"
    for path in "${paths[@]}"; do
        echo "$path"
        git -C "$path" push --force
    done
}

function squash_commits() {
    local message="$1"
    local paths=("${@:2:$#-2}")
    local time_window="${@: -1}"
    # If no message is provided, use the current date as the default value
    if [[ -z "$message" ]]; then
        message=$(date)
    fi
    # If no third argument is provided, use "1 day ago" as the default value
    if [ -n "$time_window" ]; then
        time_window="1 day ago"
    fi
    special_print "Squash"	
	for path in "${paths[@]}"; do
		latest_commit_hash=$(git -C "$path" rev-parse HEAD)
		after_time_hash=$(git -C "$path" rev-list --after="$time_window" --reverse HEAD | head -n 1)
		before_time_hash=$(git -C "$path" rev-list -n 1 --before="$time_window" HEAD)	
		if [ "$after_time_hash" != "$latest_commit_hash" ]; then
			echo "$path"
			git -C "$path" reset --soft "$before_time_hash"
			git -C "$path" commit -m "$message"
		fi
	done
	latest_commit_hash=$(git rev-parse HEAD)
	after_time_hash=$(git rev-list --after="$time_window" --reverse HEAD | head -n 1)
	before_time_hash=$(git rev-list -n 1 --before="$time_window" HEAD)
	if [ "$after_time_hash" != "$latest_commit_hash" ]; then
		echo "Parent directory"
		git reset --soft "$before_time_hash"
		git commit -m "$message"
	fi
}

function commit_paths() {
    local message="$1"
    local paths=("${@:2}")
	special_print "Commit"
    for path in "${paths[@]}"; do
		# Check if there are changes to commit
        if ! git -C "$path" diff --exit-code --quiet; then
			echo "$path"
			git -C "$path" add .
            git -C "$path" commit -m "$message"
        fi
    done
	# Check if there are changes to commit in the parent module
    if ! git diff --exit-code --quiet; then
		echo "Parent directory"
		git add .
        git commit -m "$message"
    fi
}

# Looping through all the possible entries
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help|help )
      usage
      ;;
    -ch|--checkout|checkout) 
      if [[ -n "$2" ]]; then
		checkout_paths "$2" "${paths[@]}"
        shift
      else
        echo "Error: rebase requires a valid branch name"
        exit 1
      fi
      ;;
	-pl|--pull|pull) 
		pull_paths "${paths[@]}"
		shift
      ;;
	-ps|--push|push) 
		push_paths "${paths[@]}"
		shift
      ;;
	-fp|--force-push|force-push) 
		force_push_paths "${paths[@]}"
		shift
      ;; 
	-c|--commit|commit)
		if [[ -n "$2" ]]; then
		commit_paths "$2" "${paths[@]}"
        shift
      else
        echo "Error: commit requires a comment"
        exit 1
      fi
      ;;
	-s|--squash|squash)
      if [[ -n "$2" ]]; then
		if [[ -n "$3" ]]; then
			squash_commits "$2" "${paths[@]}" "$3"
			shift
		else
			squash_commits "$2" "${paths[@]}"
			shift
		fi
      else
        echo "Error: squash requires a comment"
        exit 1
      fi
      ;;
    *)
      usage
      exit 1
      ;;
  esac
  shift
done