# Bash snippets
***
See the [bash cheatsheet gist](https://gist.github.com/LeCoupa/122b12050f5fb267e75f) for reference.
***

Logging in a pipe

```bash
# Logging function
function _log() {
	local program
	program=${0##*/}
	local now
	now=$(date '+%Y-%m-%d %H:%M:%S.%3N')

	# If stdin is not from the terminal, consume it.
	# This prevents issues with back-to-back tees.
	# It also facilitates logging outside of a pipe.
	[[ -t 0 ]] || cat > /dev/null

	echo -e "[${now}] (${program})" "$@"

} >&2
```
```bash
# Usage
source path/to/script/containing/_log.sh

# Quote "$(<command>)" for best effect
# Upstream commands ...
| tee >(_log "<message>" "$(<command>)") \
# ... downstream commands
```

Prefix a TSV with column numbers
```bash
function col_number(){
	awk 'NR==1 { for (i=1;i<=NF;i++) printf "%d%s", i, (i<NF?OFS:ORS) } 1'
}
```

Timestamps

```bash
date +"%y-%m-%d %T" # My preference
date -Iseconds # ISO 8601 (seconds precision)
```

Tidy up a temporary directory on exit

```bash
DIR_TMP=$(mktemp -d -p data)
trap 'rm -rf -- "${DIR_TMP}"' EXIT
```

Do something on exit... But not if the script completed succesfully.

```bash
trap '[[ $? -ne 0 ]] && do something' EXIT
```

Keep only filenames with find

```bash
# Only with gnu find
find <dir> -printf "%f\n"
```

Remove blank lines

```
awk 'NF'
```

Find files last modified more than 60 minutes ago

```
find <PATH> -mmin +60
```

From the `find` documentation (`man find`):

```
A numeric argument n can be specified to tests (like -amin, -mtime,  -gid,  -inum,
-links,  -size,  -uid  and -used) as
	+n     for greater than n,
	-n     for less than n,
	 n     for exactly n.
```

## Gotchas

### Grouping commands

Grouping commands in parentheses runs them in a subshell.
**Variable assignments do not persist once the subshell completes.**

```bash
( VAR=1 ); [[ -n $VAR ]]; echo $? # 1; $VAR is unset.
```

Grouping commands in braces runs them in the current shell context.
**The closing braces must be separated from the commands by a semicolon and whitespace.**

```bash
{ echo "test" }  # Fails
{ echo "test"; } # Succeeds
```

### Expressions
Prefer ```[[ expression ]]``` for text-based expressions.\
Prefer ```(( expression ))``` for arithmetic expressions.

**```$``` signs are not required for variables within these expressions.**

```$(( expression ))``` returns the value of the expression:

```bash
# Bad at testing expressions
$(( 4 + 3 )) # Fails with message: '7: command not found'

# Good at assignment
VAR=$(( 4 + 3 )); echo $VAR # 7
```

```(( expression ))``` returns nothing; its exit status depends on the expression (0 if non-zero, otherwise 1):

```bash
# Good at testing expressions
(( 1 > 2 )); echo $? # 1
(( 2 - 5 )); echo $? # 0
(( 2 - 2 )); echo $? # 1

# Bad at assignment
VAR=(( 4 + 3)) # Syntax error

# Although assignment is possible, it's harder to read.
(( VAR = 4 + 3 )); echo $VAR # 7 

# And fails when assigning a variable to 0.
(( VAR = 2 - 2 )); echo $VAR # Fails with message: 'attemped assignment to non-variable'
```

### Redirection
`>&2` Redirects stdout to the same place as sterr  
`&>2` Redirects stdout and stderr to a file named `2`
