----------------------
---- WINDOW RULES ----
----------------------
hl.window_rule({ match = { class = ".*" }, suppress_event = "maximize" })
hl.window_rule({ match = { class = "^(org\\.openrgb\\.OpenRGB|openrgb)$" }, workspace = "special:openrgb silent" })
hl.window_rule({ match = { class = "^(pavucontrol)$" },          float = true })
hl.window_rule({ match = { class = "^(nm-connection-editor)$" }, float = true })
hl.window_rule({ match = { title = "^(Picture-in-Picture)$" },   float = true, pin = true })

-- Fix XWayland apps
hl.window_rule({
    match = {
        class      = "^$",
        title      = "^$",
        xwayland   = true,
        float      = true,
        fullscreen = false,
        pin        = false,
    },
    no_focus = true,
})
