#!/usr/bin/env bash
# Bootstrap this mac setup on a fresh Mac.
# Safe to re-run: existing files are backed up to *.pre-dotfiles once, then symlinked.
set -euo pipefail

DOTFILES="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

link() {
  local src="$1" dst="$2"
  mkdir -p "$(dirname "$dst")"
  if [ -e "$dst" ] && [ ! -L "$dst" ]; then
    mv "$dst" "$dst.pre-dotfiles"
    echo "backed up $dst -> $dst.pre-dotfiles"
  fi
  ln -sfn "$src" "$dst"
  echo "linked $dst -> $src"
}

# --- Homebrew + apps ---
if ! command -v brew >/dev/null; then
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  eval "$(/opt/homebrew/bin/brew shellenv)"
fi
# older brew has no `trust`; newer refuses untrusted taps
brew trust nikitabobko/tap 2>/dev/null || true
brew trust dimentium/autoraise 2>/dev/null || true
brew bundle --file="$DOTFILES/Brewfile"

# --- oh-my-zsh + theme + plugins ---
if [ ! -d "$HOME/.oh-my-zsh" ]; then
  RUNZSH=no KEEP_ZSHRC=yes sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"
fi
ZSH_CUSTOM="${ZSH_CUSTOM:-$HOME/.oh-my-zsh/custom}"
[ -d "$ZSH_CUSTOM/themes/powerlevel10k" ] || \
  git clone --depth=1 https://github.com/romkatv/powerlevel10k.git "$ZSH_CUSTOM/themes/powerlevel10k"
[ -d "$ZSH_CUSTOM/plugins/zsh-autosuggestions" ] || \
  git clone --depth=1 https://github.com/zsh-users/zsh-autosuggestions.git "$ZSH_CUSTOM/plugins/zsh-autosuggestions"
[ -d "$ZSH_CUSTOM/plugins/zsh-syntax-highlighting" ] || \
  git clone --depth=1 https://github.com/zsh-users/zsh-syntax-highlighting.git "$ZSH_CUSTOM/plugins/zsh-syntax-highlighting"

# --- symlink configs ---
link "$DOTFILES/.zshrc"    "$HOME/.zshrc"
link "$DOTFILES/.zprofile" "$HOME/.zprofile"
link "$DOTFILES/.p10k.zsh" "$HOME/.p10k.zsh"
link "$DOTFILES/.tmux.conf" "$HOME/.tmux.conf"
link "$DOTFILES/.config/aerospace" "$HOME/.config/aerospace"
link "$DOTFILES/.config/ghostty"   "$HOME/.config/ghostty"

# --- AutoRaise settings (delay, focus, etc.) ---
defaults import com.sbmpost.AutoRaise "$DOTFILES/autoraise/com.sbmpost.AutoRaise.plist"

echo
echo "Done. Remaining manual steps:"
echo "  - Grant Accessibility permissions to AutoRaise and AeroSpace (System Settings > Privacy & Security)"
echo "  - Launch AutoRaise.app and AeroSpace once to start them"
echo "  - Restart your terminal (or: exec zsh)"
