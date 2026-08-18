# dot-files

home of the ricee

## Layout

- `fedora/` — main workstation (Fedora + Hyprland + caelestia-shell)
  - `.config/`, `.zshrc`, `.p10k.zsh`, … — configs, symlinked into `$HOME` by `deploy.sh`
  - `scripts/` — `startup_script.sh` (fresh-install bootstrap), `install_caelestia.sh` (build caelestia stack), `gaming_setup.sh`, `grub_setup.sh`
  - `docs/` — caelestia notes, TODO
- `orangepi/` — home server (Orange Pi 5 Plus, headless)
- `mac/` — macbook zshrc

## Usage

```bash
cd fedora && ./deploy.sh   # symlink configs into $HOME
```

Fresh install order: `scripts/startup_script.sh` (shell, editor, desktop packages, deploy) → `scripts/install_caelestia.sh` (build caelestia stack) → `scripts/gaming_setup.sh` (optional).
