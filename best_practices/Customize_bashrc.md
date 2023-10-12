<!--
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
-->

# Customize your `.bashrc` file

<!-- #region editable=true slideshow={"slide_type": ""} -->
# The `.bashrc` file
**"rc"** means **r**un **c**ommand

The `.bashrc` file is a script file that’s executed when a user logs in a Linux-like OS.
Since it is a dotted file, it is not visible via the standard `ls` command.
To view the complete listing of a directory, that command would need the `-a` (--all) option:
```
ls -a
```

<!-- #endregion -->

<!-- #region editable=true slideshow={"slide_type": ""} -->
# Other files related to `.bashrc`: `.bash_profile` and `.bash_aliases`

## Execution sequence at login:
  1. `.bash_profile` (or `.profile`)
  2. `.bashrc`
  3. `.bash_aliases`
<!-- #endregion -->

<!-- #region jp-MarkdownHeadingCollapsed=true -->
# 1. `.bash_profile`

### ~/.profile: executed by the command interpreter for login shells.
###  This file is not read by bash(1), if ~/.bash_profile or ~/.bash_login exists.

<!-- #endregion -->

```python editable=true slideshow={"slide_type": ""}
%pycat mybash_profile
```

<!-- #region editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true -->
# 2. `.bashrc`
## What should go into `.bashrc`?

  * ### General settings, new values for system constants, etc
  * ### Functions
<!-- #endregion -->

<!-- #region -->
## General settings example: command line history control:
```bash
# don't put duplicate lines or lines starting with space in the history.
# See bash(1) for more options
HISTCONTROL=ignoreboth

# append to the history file, don't overwrite it
shopt -s histappend

# for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
HISTSIZE=1000
HISTFILESIZE=2000
```

<!-- #endregion -->

## General settings example: colors
```
# enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"

    alias ls='ls --color=auto'
    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi
```

<!-- #region editable=true slideshow={"slide_type": ""} -->
## General settings example: line completion
(It's probably included.)

```bash
# enable programmable completion features (you don't need to enable
# this, if it's already enabled in /etc/bash.bashrc and /etc/profile
# sources /etc/bash.bashrc).
if ! shopt -oq posix; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    . /usr/share/bash-completion/bash_completion
  elif [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
  fi
fii
fiE=2000


<!-- #endregion -->

<!-- #region editable=true slideshow={"slide_type": ""} -->
## General settings example: prompt customization (PS1 and PS2)

```bash
if [ "$color_prompt" = yes ]; then
    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
else
    #PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
    PS1="[ ${debian_chroot:+($debian_chroot)}\u: \w ]\\$ "
    PS2="&gt; "
fi
```
<!-- #endregion -->

<!-- #region -->
## Functions examples
```bash
# these functions expect a path as argument:
}

l2(){
    ls -hFltr --color $1;
}

l2r(){
    ls -hFlr --color $1;
}

l2a(){
    ls -hFltra
```--color $1;
}
<!-- #endregion -->

<!-- #raw -->
%pycat mybashrc
<!-- #endraw -->

<!-- #region editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true -->
# 3. `.bash_aliases`

#### This is a well-named file! All your aliases should be there (to simplify your life).
#### If you choose this setup, the `.bashrc` should contain the following code, so that you can refer to aliases in your functions:
```bash
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.

if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi

# These functions, cdl & cda use aliases:
cdl(){
    cd $1; ll;
}

to(){
    cd $1; la;
}
```
Here, `ll` and `la` are one of the many aliases defined in my `.bash_aliases` file
<!-- #endregion -->

```python editable=true slideshow={"slide_type": ""} jupyter={"outputs_hidden": true}
%pycat mybash_aliases
```

<!-- #region editable=true slideshow={"slide_type": ""} -->
### Some options for `ls`:

__-F__: adds symbols to the end of each entry to indicate the type of the file. Here's what it does:

Files: Regular files are displayed without any additional symbols.
Directories: Directories are marked with a trailing /.
Executable Files: Executable files are marked with a trailing *.
Symbolic Links: Symbolic links are marked with a trailing @.
Sockets: Sockets are marked with a trailing =.
FIFOs (Named Pipes): FIFOs are marked with a trailing |.

__-h__: makes the file sizes in the output more human-readable, using units like KB, MB, GB, etc.
<!-- #endregion -->

<!-- #region editable=true slideshow={"slide_type": ""} -->
### Remember to source or "dot" `.bashrc` whenever you amend it or `.bash_aliases`:
```
. .bashrc
```
OR:
```
source .bashrc
```
<!-- #endregion -->
