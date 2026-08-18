# TODO

## Login screen (GDM/SDDM theming)

Currently using GDM (Fedora's default) as the display manager — it's clean but unstyled.

**Option A — Theme SDDM** (switch back to SDDM with a proper theme)

Popular themes:
- [Sugar Candy](https://github.com/Kangie/sddm-sugar-candy) — polished, highly customisable, supports background image
- [Astronaut](https://github.com/Keyitdev/sddm-astronaut-theme) — anime-styled, matches the caelestia aesthetic well
- [Where Is My SDDM Theme](https://github.com/stepanzubkov/where-is-my-sddm-theme) — minimal, blurred background

To switch back to SDDM:
```bash
sudo systemctl disable gdm
sudo systemctl enable sddm
```

Then drop the theme folder into `/usr/share/sddm/themes/` and set it in `/etc/sddm.conf`:
```ini
[Theme]
Current=<theme-name>
```

**Option B — Theme GDM**

GDM theming is limited (GNOME intentionally restricts it), but you can:
- Change the background via GNOME Shell CSS override (hacky, breaks on updates)
- Use `gdm-settings` GUI tool: `sudo dnf install gdm-settings`

**Notes**
- The caelestia lock screen (triggered by `caelestia lock` / `Super+Escape`) is separate from the login screen and already looks great.
- UWSM + Hyprland works with both GDM and SDDM — no session breakage either way.
- Whichever DM you use, logout is handled by `uwsm stop` (already set in `shell.json`).
