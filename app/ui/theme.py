"""Sistema de disenio OMNIPLAST: tipografia Sora/Inter/JetBrains Mono y paleta oceano."""

CSS = """
@import url('https://fonts.googleapis.com/css2?family=Sora:wght@600;700;800&family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap');

:root {
  --omp-navy: #0A2540; --omp-navy2: #123A5C; --omp-teal: #00A896;
  --omp-aqua: #05BFDB; --omp-coral: #FF6B6B; --omp-amber: #E9B44C;
  --omp-green: #3E9B4F; --omp-ink: #0E2433; --omp-muted: #5B7186;
  --omp-line: #E1E8F0; --omp-sand: #F5F8FB;
}

html, body, [class*="css"], .stApp, button, input, select, textarea {
  font-family: 'Inter', 'Segoe UI', system-ui, sans-serif;
}
h1, h2, h3, .omp-header h1 { font-family: 'Sora', 'Inter', sans-serif !important; letter-spacing: -0.015em; }
code, pre, [data-testid="stMetricValue"], .omp-mono { font-family: 'JetBrains Mono', Consolas, monospace; }

#MainMenu, footer { visibility: hidden; }
header[data-testid="stHeader"] { background: transparent; }

[data-testid="stMainBlockContainer"] { max-width: 1180px; padding-top: 1.4rem; }
section[data-testid="stSidebar"] {
  background: linear-gradient(180deg, #0A2540 0%, #0E2F4E 55%, #0B3B54 100%);
  border-right: 1px solid rgba(255,255,255,.06);
}
section[data-testid="stSidebar"] * { color: #DCE9F2 !important; }
section[data-testid="stSidebar"] .stRadio label, section[data-testid="stSidebar"] p {
  color: #B9CFDE !important; font-size: .88rem;
}
section[data-testid="stSidebar"] hr { border-color: rgba(255,255,255,.12); }

.stApp { background: var(--omp-sand); }

.omp-header { border-left: 5px solid var(--omp-teal); padding: .2rem 0 .2rem 1rem; margin: .1rem 0 1.1rem 0; }
.omp-header h1 { margin: 0; color: var(--omp-navy); font-size: 1.72rem; font-weight: 700; }
.omp-header .sub { color: var(--omp-muted); font-size: .93rem; margin-top: .15rem; }

.omp-hero {
  background: linear-gradient(135deg, #081E36 0%, #123A5C 55%, #0B6F7D 100%);
  border-radius: 18px; padding: 2.1rem 2.4rem; color: #EAF6FF; margin-bottom: 1.4rem;
  position: relative; overflow: hidden;
}
.omp-hero h1 { color: #FFFFFF; font-size: 2.15rem; margin: 0 0 .4rem 0; font-weight: 800; }
.omp-hero .accent { color: #4CE0C4; }
.omp-hero p { color: #BBD8E8; font-size: 1.02rem; max-width: 760px; }
.omp-chip {
  display: inline-block; background: rgba(76,224,196,.14); color: #4CE0C4;
  border: 1px solid rgba(76,224,196,.35); border-radius: 999px;
  padding: .18rem .8rem; font-size: .74rem; font-weight: 600; margin: .15rem .25rem .15rem 0;
}

.omp-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(215px, 1fr)); gap: .85rem; margin: .6rem 0 1.1rem 0; }
.omp-card {
  background: #FFFFFF; border: 1px solid var(--omp-line); border-radius: 14px;
  padding: 1.05rem 1.2rem; box-shadow: 0 1px 3px rgba(10,37,64,.05);
}
.omp-card .kpi-label { color: var(--omp-muted); font-size: .78rem; font-weight: 600;
  text-transform: uppercase; letter-spacing: .06em; }
.omp-card .kpi-value { font-family: 'JetBrains Mono', monospace; font-size: 1.62rem;
  font-weight: 600; color: var(--omp-navy); line-height: 1.25; }
.omp-card .kpi-sub { color: var(--omp-muted); font-size: .8rem; margin-top: .15rem; }

.omp-card.tall { padding: 1.15rem 1.25rem; }
.omp-card h4 { margin: 0 0 .4rem 0; color: var(--omp-navy); font-family: 'Sora'; font-size: 1rem; }
.omp-card p { color: var(--omp-muted); font-size: .86rem; margin: 0; }

.omp-badge { display: inline-block; padding: 2px 11px; border-radius: 999px;
  font-size: .72rem; font-weight: 700; letter-spacing: .02em; }
.omp-badge.ok { background: #E4F5F1; color: #00785F; }
.omp-badge.warn { background: #FDF3DC; color: #8A6116; }
.omp-badge.bad { background: #FDE8E8; color: #B42323; }
.omp-badge.info { background: #E4F1FA; color: #0B5D8A; }
.omp-badge.navy { background: #E8EEF4; color: var(--omp-navy); }

.omp-agent {
  background: #fff; border: 1px solid var(--omp-line); border-left: 4px solid var(--omp-teal);
  border-radius: 12px; padding: .9rem 1.05rem; margin-bottom: .7rem;
  transition: box-shadow .15s ease;
}
.omp-agent:hover { box-shadow: 0 4px 14px rgba(10,37,64,.10); }
.omp-agent .nombre { font-weight: 700; color: var(--omp-navy); font-family: 'Sora'; font-size: .98rem; }
.omp-agent .desc { color: var(--omp-muted); font-size: .83rem; margin: .18rem 0 .4rem 0; }

.stButton > button {
  background: var(--omp-teal); color: #fff; border: none; border-radius: 9px;
  font-weight: 600; padding: .38rem 1rem; transition: background .15s ease;
}
.stButton > button:hover { background: #008A7C; color: #fff; }
.stButton > button[kind="secondary"] { background: #EDF2F7; color: var(--omp-navy); border: 1px solid var(--omp-line); }
.stButton > button[kind="secondary"]:hover { background: #E1E9F1; }

[data-testid="stMetric"] {
  background: #fff; border: 1px solid var(--omp-line); border-radius: 12px;
  padding: .75rem .9rem;
}
[data-testid="stMetricValue"] { font-family: 'JetBrains Mono', monospace; color: var(--omp-navy); }

[data-testid="stDataFrame"] { border: 1px solid var(--omp-line); border-radius: 10px; }

.stTabs [data-baseweb="tab-list"] { gap: 2px; border-bottom: 2px solid var(--omp-line); }
.stTabs [data-baseweb="tab"] {
  border-radius: 8px 8px 0 0; padding: .45rem 1rem; font-weight: 600; color: var(--omp-muted);
  background: transparent;
}
.stTabs [aria-selected="true"] { background: #fff; color: var(--omp-navy); }
.stTabs [data-baseweb="tab-highlight"] { background-color: var(--omp-teal) !important; height: 3px; }

div[data-testid="stExpander"] { border: 1px solid var(--omp-line); border-radius: 12px; }

.omp-footer { margin-top: 2.2rem; padding-top: 1rem; border-top: 1px solid var(--omp-line);
  color: var(--omp-muted); font-size: .78rem; }
.omp-footer a { color: var(--omp-teal); text-decoration: none; }

.omp-section-title { display:flex; align-items:center; gap:.6rem; margin: 1.6rem 0 .8rem 0;
  color: var(--omp-navy); font-family: 'Sora'; font-size: 1.12rem; font-weight: 700; }
.omp-section-title .bar { width: 26px; height: 5px; border-radius: 3px; background: var(--omp-teal); display: inline-block; }

.omp-resumen {
  background: linear-gradient(90deg, #EAF7F5, #F2FAF9); border: 1px solid #BCE5DE;
  border-left: 4px solid var(--omp-teal); border-radius: 10px; padding: .85rem 1rem;
  color: #073B34; font-size: .92rem;
}
.omp-requiere {
  background: #FDF6E7; border: 1px solid #F0DCB2; border-left: 4px solid var(--omp-amber);
  border-radius: 10px; padding: .85rem 1rem; color: #5C4408; font-size: .92rem;
}
.omp-error {
  background: #FDEBEB; border: 1px solid #F2C5C5; border-left: 4px solid var(--omp-coral);
  border-radius: 10px; padding: .85rem 1rem; color: #7A1F1F; font-size: .92rem;
}
.omp-refs { font-size: .78rem; color: var(--omp-muted); }
.omp-refs span { display: block; margin-bottom: .12rem; }
"""
