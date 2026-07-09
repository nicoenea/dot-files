# mac

macOS setup: AeroSpace + AutoRaise + Ghostty (with shaders) + oh-my-zsh/powerlevel10k.

## What's here

| Path | What it is |
|---|---|
| `.zshrc`, `.zprofile`, `.p10k.zsh` | zsh config |
| `.config/aerospace/` | AeroSpace tiling WM config |
| `.config/ghostty/` | Ghostty config + custom shaders (balatro, crt) |
| `autoraise/` | AutoRaise settings exported from macOS defaults |
| `Brewfile` | Apps and fonts |
| `install.sh` | One-shot bootstrap for a fresh Mac |

## Setup on a new Mac

```sh
git clone https://github.com/nicoenea/dot-files.git ~/repos/dot-files
cd ~/repos/dot-files/mac && ./install.sh
```

The script installs Homebrew (if missing), the Brewfile apps, oh-my-zsh + powerlevel10k + plugins, symlinks all configs into place, and imports the AutoRaise settings. It backs up any existing files to `*.pre-dotfiles` and is safe to re-run.

Manual follow-ups after install: grant Accessibility permissions to AutoRaise and AeroSpace, launch each once, restart the terminal.

## AutoRaise settings

AutoRaise settings live in macOS defaults, not in a file: after changing them, re-export with

```sh
defaults export com.sbmpost.AutoRaise autoraise/com.sbmpost.AutoRaise.plist
```
