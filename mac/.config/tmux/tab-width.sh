#!/bin/sh
# Resize tmux "tabs" to share the status bar equally, like native macOS tabs.
# Called from hooks in tmux.conf on attach/resize/tab-open/tab-close.
S="$1"
[ -n "$S" ] || S="$(tmux display -p '#{session_id}')" || exit 0

W="$(tmux display -p -t "$S" '#{client_width}')"
N="$(tmux display -p -t "$S" '#{session_windows}')"
case "$W" in ''|*[!0-9]*|0) W=120 ;; esac
case "$N" in ''|*[!0-9]*|0) exit 0 ;; esac

PAD=$(( W / N ))
[ "$PAD" -gt 60 ] && PAD=60   # don't let a lone tab span a huge monitor
[ "$PAD" -lt 12 ] && PAD=12   # keep tabs legible when there are many
INNER=$(( PAD - 3 ))          # 2 leading spaces + 1 trailing gap
TRIM=$(( INNER - 2 ))

IN="#[bg=#3b3d4a,fg=#a8a8a8]  #{p${INNER}:#{=${TRIM}:pane_title}}#[bg=default] "
CUR="#[bg=#57c7ff,fg=#282a36,bold]  #{p${INNER}:#{=${TRIM}:pane_title}}#[bg=default] "

tmux list-windows -t "$S" -F '#{window_id}' | while read -r w; do
  tmux setw -t "$w" window-status-format "$IN"
  tmux setw -t "$w" window-status-current-format "$CUR"
done
