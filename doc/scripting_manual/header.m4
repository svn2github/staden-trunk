@c ---------------------------------------------------------------------------
@c Experiment with smaller amounts of whitespace between chapters
@c and sections.
@c ---------------------------------------------------------------------------
@tex
@set tex
\global\chapheadingskip = 15pt plus 4pt minus 2pt 
\global\secheadingskip = 12pt plus 3pt minus 2pt
\global\subsecheadingskip = 9pt plus 2pt minus 2pt
@end tex

@c ---------------------------------------------------------------------------
@c @split{} command
@c
@c only makes sense for html.
@c ---------------------------------------------------------------------------
@tex
\global\def\split{}
@end tex

@c ---------------------------------------------------------------------------
@c Experiment with smaller amounts of whitespace between paragraphs in
@c the 8.5 by 11 inch `format'.
@tex
\global\parskip 6pt plus 1pt
@end tex
@c ---------------------------------------------------------------------------

@c ---------------------------------------------------------------------------
@c Magic with comments. m4 can set comment characters to whatever it wants.
@c They do not even have to be on one line (but by default the start and end
@c characters are "#" and newline).
@c
@c We `define' new start and end comments: @nm4 and @m4. (Remember as no m4 and
@c m4).
@c
@c m4 will not remove text in comments, it just ignores it. So the comment
@c characters themselves need to be harmless to tex. We solve this by creating
@c two new tex commands to do nothing.
@c ---------------------------------------------------------------------------
@tex
\global\def\m4{}
\global\def\nm4{}
@end tex
changecom(@nm4,@m4)

@c ---------------------------------------------------------------------------
@c Rename the m4 commands to _commands. This will greatly reduce the chance of
@c them occurring in our text by chance.
@c ---------------------------------------------------------------------------
define(`_define',defn(`define'))
define(`_changecom',defn(`changecom'))
define(`_changequote',defn(`changequote'))
define(`_errprint',defn(`errprint'))
define(`_maketemp',defn(`maketemp'))
define(`_sinclude',defn(`sinclude'))
define(`_translit',defn(`translit'))
define(`_traceoff',defn(`traceoff'))
define(`_undefine',defn(`undefine'))
define(`_undivert',defn(`undivert'))
define(`_decr',defn(`decr'))
define(`_defn',defn(`defn'))
define(`_divert',defn(`divert'))
define(`_divnum',defn(`divnum'))
define(`_dlen',defn(`dlen'))
define(`_dumpdef',defn(`dumpdef'))
define(`_eval',defn(`eval'))
define(`_m4exit',defn(`m4exit'))
define(`_ifelse',defn(`ifelse'))
define(`_ifdef',defn(`ifdef'))
define(`_include',defn(`include'))
define(`_incr',defn(`incr'))
define(`_index',defn(`index'))
define(`_popdef',defn(`popdef'))
define(`_pushdef',defn(`pushdef'))
define(`_shift',defn(`shift'))
define(`_substr',defn(`substr'))
define(`_syscmd',defn(`syscmd'))
define(`_sysval',defn(`sysval'))
define(`_traceon',defn(`traceon'))
define(`_m4wrap',defn(`m4wrap'))
define(`_format',define(`format'))

_undefine(`define')
_undefine(`changecom')
_undefine(`changequote')
_undefine(`errprint')
_undefine(`maketemp')
_undefine(`sinclude')
_undefine(`translit')
_undefine(`traceoff')
_undefine(`undefine')
_undefine(`undivert')
_undefine(`unix')
_undefine(`windows')
_undefine(`decr')
_undefine(`defn')
_undefine(`divert')
_undefine(`divnum')
_undefine(`dlen')
_undefine(`dumpdef')
_undefine(`eval')
_undefine(`m4exit')
_undefine(`ifelse')
_undefine(`ifdef')
_undefine(`include')
_undefine(`incr')
_undefine(`index')
_undefine(`popdef')
_undefine(`pushdef')
_undefine(`shift')
_undefine(`substr')
_undefine(`syscmd')
_undefine(`sysval')
_undefine(`traceon')
_undefine(`m4wrap')
_undefine(`format')

@c ---------------------------------------------------------------------------
@c Change quotes to [[ and ]]. Otherwise quotes are likely to cause us
@c problems. [[ and ]] are not likely to occur by chance in our docs.
@c
@c If we need to use an m4 keyword in our text, then we may do so with
@c (eg) [[_m4command]].
@c
@c If we wish to use [[ and ]] in our text, enclose it with comments:
@c @nm4{}[[@m4{}
@c ---------------------------------------------------------------------------
_changequote([[,]])

@c ---------------------------------------------------------------------------
@c picture macro
@c
@c Adds a picture to the document. For tex it loads a PostScript file. For
@c html it loads a gif file.
@c
@c argument 1: a filename prefix. .ps and .gif are added to the prefix
@c             as required.
@c ---------------------------------------------------------------------------
_define([[_picture]],[[_ifdef([[_unix]],[[_ifdef([[_tex]],[[@tex
@sp 1
@epsfbox{[[$*]].unix.ps}
@end tex]])
_ifdef([[_html]],[[
@ifhtml
<p>
<img src="[[$*]].unix.gif" alt="[picture]">
@end ifhtml]])]],[[_ifdef([[_tex]],[[@tex
@sp 1
@epsfbox{[[$*]].unix.ps}
@end tex]])
_ifdef([[_html]],[[
@ifhtml
<p>
<img src="[[$*]].unix.gif" alt="[picture]">
@end ifhtml]])]])]])

@c ---------------------------------------------------------------------------
@c lpicture macro
@c
@c Adds a large picture to the document. In tex this is the same as the
@c picture macro. For html it displays a small gif file with a link to the
@c full size one.
@c
@c argument 1: a filename prefix. .ps, .gif, .small.gif and .gif.html are
@c             added to the prefix as required.
@c ---------------------------------------------------------------------------
_define([[_lpicture]],[[_ifdef([[_unix]],[[_ifdef([[_tex]],[[@tex
@sp 1
@epsfbox{[[$*]].unix.ps}
@end tex]])
_ifdef([[_html]],[[
@ifhtml
<p>
<a href="[[$*]].unix.gif.html"><img src="[[$*]].small.unix.gif" alt="[picture]"></a>
<br><font size="-1">(Click for full size image)<font size="+0"><br>
@end ifhtml]])]],[[_ifdef([[_tex]],[[@tex
@sp 1
@epsfbox{[[$*]].unix.ps}
@end tex]])
_ifdef([[_html]],[[
@ifhtml
<p>
<a href="[[$*]].unix.gif.html"><img src="[[$*]].small.unix.gif" alt="[picture]"></a>
<br><font size="-1">(Click for full size image)<font size="+0"><br>
@end ifhtml]])]])]])

@c ---------------------------------------------------------------------------
@c @nm4{}
@c
@c _ifunix macro
@c _ifwindows macro
@c
@c These two macros may be used to surround text which we wish to only
@c appear in one version or another. They check the _ifunix and _ifwindows
@c defines.
@c An example usage is:
@c
@c     _ifunix([[
@c     @split{}
@c     @node Assembly-CAP2
@c     @section Assembly CAP2
@c     _include(cap2-t.texi)
@c     ]])(
@c
@c An alternative to this is using _ifdef directly. Eg:
@c
@c     _ifdef([[_unix]],[[
@c     @split{}
@c     @node Assembly-CAP2
@c     @section Assembly CAP2
@c     _include(cap2-t.texi)
@c     ]])(
@c
@c @m4{}
@c ---------------------------------------------------------------------------
_define([[_ifunix]],[[_ifdef([[_unix]],[[$*]])]])
_define([[_ifwindows]],[[_ifdef([[_windows]],[[$*]])]])

@c ---------------------------------------------------------------------------_
@c uref macro
@c
@c This exists in newer texinfo release, but for now we try to emulate it as
@c well as possible (albeit in a m4 instead of texinfo manner).
@c
@c _uref(url) will just link to that url, with the 'url' as the text.
@c _uref(url,text) will link to that url, with 'text' as the text in the
@c   html format. For tex format it'll use "text (@code{url})".
@c _uref(url,,text) will link to that url, with 'text' as the text in both
@c   html and tex formats.
@c ---------------------------------------------------------------------------
_define([[_uref]],[[_ifelse(1,$#,[[_ifdef([[_html]],[[@ifhtml
<a href="$1">
@end ifhtml]])
$1
_ifdef([[_html]],[[
@ifhtml
</a>
@end ifhtml]])]],[[_ifelse(2,$#,[[_ifdef([[_tex]],[[$2 (@code{$1})]])
_ifdef([[_html]],[[@ifhtml
<a href="$1">$2</a>
@end ifhtml]])]],[[_ifdef([[_html]],[[
@ifhtml

<a href="$1">
@end ifhtml
$3
@ifhtml
</a>
@end ifhtml
]])]])]])]])

@c normal refs
_ifdef([[_tex]],[[
_define([[_fxref]],[[@xref{$1,$1,$2}.]])
_define([[_fpref]],[[@pxref{$1,$1,$2}]])
_define([[_fref]],[[@ref{$1,$1,$2}.]])
_define([[_split]],[[]])
]])

@c html refs
_ifdef([[_html]],[[
_define([[_fxref]],[[
@ifhtml
<!-- XREF:$1 -->
@end ifhtml
@xref{$1,$1,$2,$3,$3}.]])
_define([[_fpref]],[[
@ifhtml
<!-- XREF:$1 -->
@end ifhtml
@pxref{$1,$1,$2,$3,$3}]])
_define([[_fref]],[[
@ifhtml
<!-- XREF:$1 -->
@end ifhtml
@ref{$1,$1,$2,$3,$3}.]])
_define([[_split]],[[@split]])
]])

@c common refs
_define([[_oxref]],[[@xref{$1,$1,$2}]])
_define([[_oref]],[[@ref{$1,$1,$2}]])

@c A horizontal ruler, using TeX or HTML
_define([[_rule]],[[@sp 1
@tex
\hrule height 0.5pt width \hsize
@end tex
@ifhtml
<hr>
@end ifhtml]])
