# System Configuration — Nic's Fedora Workstation

## OS & Environment
- **Distro:** Fedora 43
- **Desktop:** Hyprland on **Wayland**
- **GPU:** NVIDIA RTX 3080 (proprietary driver; NVIDIA env vars set in `hyprland.conf`)
- **Session type:** Wayland — do not suggest X11-only solutions
- **Package manager:** `dnf`

## Hardware
- **Monitors:**
  - `DP-2` — ASUS VG27AC1A, left/main, 2560x1440 @ 170Hz — workspaces 1–5
  - `DP-3` — Gigabyte M27Q, right, 2560x1440 @ 170Hz — workspaces 6–10
- **Resolution:** 2K (2560x1440)

## Desktop Environment

### Hyprland config
**File:** `.config/hypr/hyprland.conf`
- Layout: dwindle
- Mod key: `Super`
- Terminal: `ghostty`
- Launcher: `caelestia shell drawers toggle launcher`

### caelestia-shell
Replaces both Waybar (bar) and Mako (notifications). Also provides: OSD overlays, app launcher, lock screen, session menu.
- Autostarted via `exec-once = caelestia shell -d` in hyprland.conf
- Config: `~/.config/caelestia/shell.json` (symlinked from dotfiles)
- Bar excluded from DP-3 (right monitor) via `shell.json`
- See `caelestia.md` for full details on the caelestia setup

> **Note:** Waybar and Mako configs still exist in this repo but are NOT active.
> They're kept for reference only (deploy.sh has their symlinks commented out).

### Wallpaper
`swww-daemon` runs at startup. Use `caelestia wallpaper set <file>` to change wallpaper.
Wallpapers live in `~/Pictures/Wallpapers`.

### Clipboard
`cliphist` manages clipboard history. Access via `Super+;`.

## Keybindings (Hyprland)

### Core
| Key | Action |
|-----|--------|
| `Super+Return` | Open terminal (ghostty) |
| `Super+Space` | Toggle app launcher |
| `Super+C` | Kill active window |
| `Super+F` | Fullscreen |
| `Super+Shift+F` | Fullscreen (maximize, keep bar) |
| `Super+V` | Toggle floating |
| `Super+G` | Game mode (suppresses notifications) |
| `Super+Escape` | Lock screen |
| `Super+Shift+Ctrl+Q` | Exit Hyprland |

### Focus & Window Management
| Key | Action |
|-----|--------|
| `Super+H/J/K/L` | Move focus left/down/up/right |
| `Super+Shift+H/J/K/L` | Move window left/down/up/right |
| `Super+Ctrl+Arrows` | Resize active window (40px steps) |
| `Super+,` / `Super+.` | Focus left/right monitor |
| `Super+Shift+,` / `Super+Shift+.` | Move workspace to left/right monitor |

### Workspaces
| Key | Action |
|-----|--------|
| `Alt+1-5` | Switch to workspace 1–5 (left monitor, DP-2) |
| `Alt+Q/W/E/R/T` | Switch to workspace 6–10 (right monitor, DP-3) |
| `Shift+Alt+1-5` | Move window to workspace 1–5 |
| `Shift+Alt+Q/W/E/R/T` | Move window to workspace 6–10 |
| `Super+Scroll` | Cycle through workspaces |

### Mouse
| Action | Binding |
|--------|---------|
| Move window | `Super+LMB drag` |
| Resize window | `Super+RMB drag` |

### Utilities
| Key | Action |
|-----|--------|
| `Ctrl+Shift+3` | Screenshot full screen → clipboard |
| `Ctrl+Shift+4` | Screenshot region → clipboard |
| `Super+;` | Clipboard history (cliphist + wofi) |
| Media keys | Volume (wpctl), playback (playerctl), brightness (brightnessctl) |

## Shell & Terminal
- **Shell:** zsh with oh-my-zsh + Powerlevel10k theme
- **Terminal:** Ghostty (custom shaders: balatro.glsl, crt.glsl)
- **Prompt:** Powerlevel10k (`~/.p10k.zsh`)
- **Notable tools:** `eza` (modern ls, via cargo), `starship` (via install script), `btop`, `fastfetch`

## Dotfile Management
**Script:** `deploy.sh` — symlinks all configs into `$HOME`. Run as regular user.
**Install caelestia:** `install_caelestia.sh` — builds quickshell + caelestia-shell from source on Fedora.

## Notes for Claude Code
- This is a **Wayland** system — always prefer Wayland-compatible solutions
- **NVIDIA + Wayland**: env vars like `LIBVA_DRIVER_NAME=nvidia`, `__GLX_VENDOR_LIBRARY_NAME=nvidia`, `NVD_BACKEND=direct` are set in hyprland.conf
- Package manager is `dnf` — confirm package names exist in Fedora repos before recommending
- caelestia-shell is the bar/notification daemon — don't suggest waybar or mako as replacements
- The `caelestia.md` file in this repo has detailed notes on the caelestia install, architecture, and known Fedora-specific issues
