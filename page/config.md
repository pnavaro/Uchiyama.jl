@def title       = "Uchiyama.jl"
@def prepath     = "Uchiyama.jl"
@def description = """
                   Julia implementation of the Uchiyama particle model
                   """
@def authors     = "Nathalie Ayi and Pierre Navaro"

<!--  NAVBAR SPECS
  NOTE:
  - add_docs:  whether to add a pointer to your docs website
  - docs_url:  the url of the docs website (ignored if add_docs=false)
  - docs_name: how the link should be named in the navbar

  - add_nav_logo:  whether to add a logo left of the package name
  - nav_logo_path: where the logo is
-->
@def add_docs  = false
@def docs_url  = "https://franklinjl.org/"
@def docs_name = "Docs"

@def add_nav_logo   = true
@def nav_logo_path  = "/assets/logo.svg"
@def nav_logo_alt   = "Logo"
@def nav_logo_style = """
                      height:         25px;
                      padding-right:  10px;
                      """

@def use_header         = true
@def use_header_img     = true
@def header_img_path    = "url(\"assets/diagonal-lines.svg\")"
@def header_img_style   = """
                          background-repeat: repeat;
                          """
@def header_margin_top  = "55px" <!-- 55-60px ~ touching nav bar -->

@def use_hero           = false
@def hero_width         = "80%"
@def hero_margin_top    = "100px"

@def add_github_view  = true
@def add_github_star  = true
@def github_repo      = "pnavaro/Uchiyama.jl"

<!-- SECTION LAYOUT
NOTE:
  - section_width:  integer number to control the default width of sections
                    you can also set it for individual sections by specifying
                    the width argument: `\begin{:section, ..., width=10}`.
-->
@def section_width = 10

<!-- COLOR PALETTE
You can use Hex, RGB or SVG color names; these tools are useful to choose:
  - color wheel: https://developer.mozilla.org/en-US/docs/Web/CSS/CSS_Colors/Color_picker_tool
  - color names: https://developer.mozilla.org/en-US/docs/Web/CSS/color_value

NOTE:
  - header_color:      background color of the header
  - link_color:        color of links
  - link_hover_color:  color of links when hovered
  - section_bg_color:  background color of "secondary" sections to help
                       visually separate between sections.
  - footer_link_color: color of links in the footer
-->
@def header_color       = "#3f6388"
@def link_color         = "#2669DD"
@def link_hover_color   = "teal"
@def section_bg_color   = "#f6f8fa"
@def footer_link_color  = "cornflowerblue"

<!-- CODE LAYOUT
NOTE:
  - highlight_theme:    theme for the code, pick one from
                        https://highlightjs.org/static/demo/ for instance
                        "github" or "atom-one-dark"; use lower case and replace
                        spaces with `-`.
  - code_border_radius: how rounded the corners of code blocks should be
  - code_output_indent: how much left-identation to add for "output blocks"
                        (results of the evaluation of code blocks), use 0 if
                        you don't want indentation.
-->
@def highlight_theme    = "atom-one-dark"
@def code_border_radius = "10px"
@def code_output_indent = "15px"


<!-- YOUR DEFINITIONS
See franklinjl.org for more information on how to introduce your own
definitions and how they can be useful.
-->


<!-- INTERNAL DEFINITIONS =====================================================
===============================================================================
These definitions are important for the good functioning of some of the
commands that are defined and used in PkgPage.jl
-->
@def sections        = Pair{String,String}[]
@def section_counter = 1
@def showall         = true

\newcommand{\html}[1]{~~~#1~~~}

\newenvironment{center}{\html{<div style="text-align:center;">}}{\html{</div>}}

\newenvironment{columns}{\html{<div class="container"><div class="row">}}{\html{</div></div>}}
