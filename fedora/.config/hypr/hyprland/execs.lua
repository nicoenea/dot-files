-------------------
---- AUTOSTART ----
-------------------
hl.on("hyprland.start", function()
    hl.exec_cmd("gnome-keyring-daemon --start --components=secrets,ssh,gpg")
    hl.exec_cmd("swww-daemon")
    -- RGB lighting: HyperHDR drives the Skydimo strip (native Skydimo mode, EGL/DMA screen capture)
    -- OpenRGB handles all other devices (RAM, GPU, Kraken, mouse, board) — no longer controls the strip
    hl.exec_cmd("bash -c 'sleep 3; skydimo-init && hyperhdr'")
    hl.exec_cmd("bash -c 'sleep 3; openrgb --gui'")
    -- caelestia-shell replaces both waybar (bar) and mako (notifications).
    -- It also handles: OSD overlays, app launcher, lock screen, session menu.
    -- Wait for swww socket before starting caelestia so wallpaper loads consistently.
    hl.exec_cmd("bash -c 'until swww query &>/dev/null; do sleep 0.1; done; swww restore; caelestia shell -d'")
    hl.exec_cmd("sh -c '/usr/lib/polkit-gnome/polkit-gnome-authentication-agent-1 || /usr/libexec/polkit-gnome-authentication-agent-1'")
    hl.exec_cmd("nm-applet --indicator")
    hl.exec_cmd("wl-paste --type text --watch cliphist store")
    hl.exec_cmd("wl-paste --type image --watch cliphist store")
    hl.exec_cmd("ydotoold")
end)
