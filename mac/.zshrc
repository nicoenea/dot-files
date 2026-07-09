# Ghostty runs its tabs as tmux windows: macOS-native tabs register as
# separate windows with tiling WMs (AeroSpace), breaking layouts.
if [ -z "$TMUX" ] && [ "$TERM_PROGRAM" = "ghostty" ] && [[ -o interactive ]] \
   && command -v tmux >/dev/null; then
  exec tmux new-session
fi

# Enable Powerlevel10k instant prompt. Should stay close to the top of ~/.zshrc.
# Initialization code that may require console input (password prompts, [y/n]
# confirmations, etc.) must go above this block; everything else may go below.
if [[ -r "${XDG_CACHE_HOME:-$HOME/.cache}/p10k-instant-prompt-${(%):-%n}.zsh" ]]; then
  source "${XDG_CACHE_HOME:-$HOME/.cache}/p10k-instant-prompt-${(%):-%n}.zsh"
fi

# If you come from bash you might have to change your $PATH.
# export PATH=$HOME/bin:/usr/local/bin:$PATH

# Path to your oh-my-zsh installation.
export ZSH="$HOME/.oh-my-zsh"

# Set name of the theme to load --- if set to "random", it will
# load a random theme each time oh-my-zsh is loaded, in which case,
# to know which specific one was loaded, run: echo $RANDOM_THEME
# See https://github.com/ohmyzsh/ohmyzsh/wiki/Themes
ZSH_THEME="powerlevel10k/powerlevel10k"
#ZSH_THEME="spaceship"

# Set list of themes to pick from when loading at random
# Setting this variable when ZSH_THEME=random will cause zsh to load
# a theme from this variable instead of looking in $ZSH/themes/
# If set to an empty array, this variable will have no effect.
# ZSH_THEME_RANDOM_CANDIDATES=( "robbyrussell" "agnoster" )

# Uncomment the following line to use case-sensitive completion.
# CASE_SENSITIVE="true"

# Uncomment the following line to use hyphen-insensitive completion.
# Case-sensitive completion must be off. _ and - will be interchangeable.
# HYPHEN_INSENSITIVE="true"

# Uncomment one of the following lines to change the auto-update behavior
# zstyle ':omz:update' mode disabled  # disable automatic updates
# zstyle ':omz:update' mode auto      # update automatically without asking
# zstyle ':omz:update' mode reminder  # just remind me to update when it's time

# Uncomment the following line to change how often to auto-update (in days).
# zstyle ':omz:update' frequency 13

# Uncomment the following line if pasting URLs and other text is messed up.
# DISABLE_MAGIC_FUNCTIONS="true"

# Uncomment the following line to disable colors in ls.
# DISABLE_LS_COLORS="true"

# Uncomment the following line to disable auto-setting terminal title.
# DISABLE_AUTO_TITLE="true"

# Uncomment the following line to enable command auto-correction.
# ENABLE_CORRECTION="true"

# Uncomment the following line to display red dots whilst waiting for completion.
# You can also set it to another string to have that shown instead of the default red dots.
# e.g. COMPLETION_WAITING_DOTS="%F{yellow}waiting...%f"
# Caution: this setting can cause issues with multiline prompts in zsh < 5.7.1 (see #5765)
# COMPLETION_WAITING_DOTS="true"

# Uncomment the following line if you want to disable marking untracked files
# under VCS as dirty. This makes repository status check for large repositories
# much, much faster.
# DISABLE_UNTRACKED_FILES_DIRTY="true"

# Uncomment the following line if you want to change the command execution time
# stamp shown in the history command output.
# You can set one of the optional three formats:
# "mm/dd/yyyy"|"dd.mm.yyyy"|"yyyy-mm-dd"
# or set a custom format using the strftime function format specifications,
# see 'man strftime' for details.
# HIST_STAMPS="mm/dd/yyyy"

# Would you like to use another custom folder than $ZSH/custom?
# ZSH_CUSTOM=./zshcustom

# Which plugins would you like to load?
# Standard plugins can be found in $ZSH/plugins/
# Custom plugins may be added to $ZSH_CUSTOM/plugins/
# Example format: plugins=(rails git textmate ruby lighthouse)
# Add wisely, as too many plugins slow down shell startup.
plugins=(
  git
  brew
  common-aliases
  node
  npm
  rand-quote
  sudo
  yarn
  z
  colored-man-pages
  colorize
  cp
  zsh-syntax-highlighting
  zsh-autosuggestions
)

source $ZSH/oh-my-zsh.sh

# .zshrc
# include Z, yo
#. ~/z.sh

# Pyenv config 
export PYENV_ROOT="$HOME/.pyenv"
command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"
command -v pyenv >/dev/null && eval "$(pyenv init -)"

# nvm config ($HOMEBREW_PREFIX comes from brew shellenv in .zprofile)
# Sourcing nvm.sh costs ~2s, so lazy-load it on first use of nvm/node/npm/npx
export NVM_DIR="$HOME/.nvm"
_load_nvm() {
  unset -f nvm node npm npx 2>/dev/null
  [ -s "$HOMEBREW_PREFIX/opt/nvm/nvm.sh" ] && \. "$HOMEBREW_PREFIX/opt/nvm/nvm.sh"
}
nvm() { _load_nvm; nvm "$@"; }
if ! command -v node >/dev/null; then
  node() { _load_nvm; node "$@"; }
  npm()  { _load_nvm; npm "$@"; }
  npx()  { _load_nvm; npx "$@"; }
fi

# asdf: pre-0.16 has asdf.sh; 0.16+ is a Go binary that only needs shims on PATH
if [ -f "$HOMEBREW_PREFIX/opt/asdf/libexec/asdf.sh" ]; then
  . "$HOMEBREW_PREFIX/opt/asdf/libexec/asdf.sh"
elif command -v asdf >/dev/null; then
  export PATH="${ASDF_DATA_DIR:-$HOME/.asdf}/shims:$PATH"
fi


# User configuration

# export MANPATH="/usr/local/man:$MANPATH"

# You may need to manually set your language environment
# export LANG=en_US.UTF-8

# Preferred editor for local and remote sessions
# if [[ -n $SSH_CONNECTION ]]; then
#   export EDITOR='vim'
# else
#   export EDITOR='mvim'
# fi

# Compilation flags
# export ARCHFLAGS="-arch x86_64"
export USE_GKE_GCLOUD_AUTH_PLUGIN=True
# Set personal aliases, overriding those provided by oh-my-zsh libs,
# plugins, and themes. Aliases can be placed here, though oh-my-zsh
# users are encouraged to define aliases within the ZSH_CUSTOM folder.
# For a full list of active aliases, run `alias`.
#
# Example aliases
alias zshconfig="code ~/.zshrc"
alias yabaiconfig="vim ~/.yabairc"
alias sketchybarrc="vim ~/.config/sketchybar/sketchybarrc"
alias skhdconfig="vim ~/.skhdrc"
alias rerice="brew services restart yabai && brew services restart sketchybar"
alias updatesettings="source ~/.zshrc"
alias yabairestart="brew services restart yabai"
alias sketchrestart="brew services restart sketchybar"
alias fastlane='bundle exec fastlane'
alias cd..='cd ..'
alias cat="bat"
# alias ohmyzsh="mate ~/.oh-my-zsh"fpath=($fpath "/Users/nicolas.enea/.zfunctions")

# To customize prompt, run `p10k configure` or edit ~/.p10k.zsh.
[[ ! -f ~/.p10k.zsh ]] || source ~/.p10k.zsh

export PATH="/Applications/Docker.app/Contents/Resources/bin:$PATH"

# The next line updates PATH for the Google Cloud SDK.
if [ -f '/Users/nicolas.enea/google-cloud-sdk/path.zsh.inc' ]; then . '/Users/nicolas.enea/google-cloud-sdk/path.zsh.inc'; fi

# The next line enables shell command completion for gcloud.
if [ -f '/Users/nicolas.enea/google-cloud-sdk/completion.zsh.inc' ]; then . '/Users/nicolas.enea/google-cloud-sdk/completion.zsh.inc'; fi
export PATH="$HOME/.local/bin:$PATH"

# Inside tmux, oh-my-zsh stops emitting xterm-style titles, so pane_title
# (what the tab bar shows) stays stuck on the hostname. Publish it ourselves:
# directory when idle, command while running — like native tabs did.
# Keep this at the very end: earlier init (p10k) rebuilds the hook arrays.
if [ -n "$TMUX" ]; then
  _tab_title_precmd()  { printf '\e]2;%s\a' "${PWD/#$HOME/~}"; }
  _tab_title_preexec() { printf '\e]2;%s\a' "$1"; }
  autoload -Uz add-zsh-hook
  add-zsh-hook precmd  _tab_title_precmd
  add-zsh-hook preexec _tab_title_preexec
fi
