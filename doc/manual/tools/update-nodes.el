(load-library "texnfo-upd")
(set-mark (point))

;; 28/02/00 jkb
;; The end-of-buffer crashes xemacs command, but the docs now tell us to
;; use goto-char instead.
;; (end-of-buffer)
(goto-char (point-max))

(texinfo-update-node)
(save-some-buffers t)
