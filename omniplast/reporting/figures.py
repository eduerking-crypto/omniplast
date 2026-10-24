"""Plantilla grafica OMNIPLAST para plotly: paleta, tipografia y layouts consistentes."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

from ..config import CHART_SEQUENCE, PALETTE

LAYOUT = dict(
    template="plotly_white",
    font=dict(family="Inter, Segoe UI, sans-serif", size=13, color=PALETTE["ink"]),
    title=dict(font=dict(family="Sora, Inter, sans-serif", size=17, color=PALETTE["navy"]),
               x=0.01, pad=dict(t=24, b=6)),
    margin=dict(l=52, r=24, t=64, b=48),
    hoverlabel=dict(font=dict(family="Inter, sans-serif", size=12)),
    paper_bgcolor="white", plot_bgcolor="white",
    colorway=CHART_SEQUENCE,
    xaxis=dict(gridcolor=PALETTE["line"], zerolinecolor=PALETTE["line"]),
    yaxis=dict(gridcolor=PALETTE["line"], zerolinecolor=PALETTE["line"]),
    legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
)


def _apply(fig, title=None, ytitle=None, xtitle=None, height=420):
    fig.update_layout(**{k: v for k, v in LAYOUT.items() if k != "template"})
    fig.update_layout(template="plotly_white", height=height)
    if title:
        fig.update_layout(title=dict(text=title, **LAYOUT["title"]))
    if ytitle:
        fig.update_yaxes(title_text=ytitle)
    if xtitle:
        fig.update_xaxes(title_text=xtitle)
    return fig


def bars(df, x, y, title=None, color=None, ytitle="", xtitle="", horizontal=False,
         height=420, text_auto=True):
    if horizontal:
        fig = px.bar(df, y=x, x=y, color=color, orientation="h",
                     color_discrete_sequence=CHART_SEQUENCE)
    else:
        fig = px.bar(df, x=x, y=y, color=color, color_discrete_sequence=CHART_SEQUENCE,
                     barmode="group")
    if text_auto:
        fig.update_traces(texttemplate="%{value:,.3s}", textposition="outside",
                          cliponaxis=False)
    return _apply(fig, title, ytitle, xtitle, height)


def pie(labels, values, title=None, height=420):
    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.42,
                                 marker=dict(colors=CHART_SEQUENCE,
                                             line=dict(color="white", width=2)),
                                 texttemplate="%{percent:.1%}")])
    return _apply(fig, title, height=height)


def box(df, x, y, title=None, color=None, ytitle="", xtitle="", height=440, points="all"):
    fig = px.box(df, x=x, y=y, color=color, points=points,
                 color_discrete_sequence=CHART_SEQUENCE)
    fig.update_traces(boxmean=True)
    return _apply(fig, title, ytitle, xtitle, height)


def scatter(df, x, y, title=None, color=None, size=None, trendline=False,
            ytitle="", xtitle="", height=440, text=None):
    fig = px.scatter(df, x=x, y=y, color=color, size=size, text=text,
                     color_discrete_sequence=CHART_SEQUENCE)
    if trendline:
        x_num = pd.to_numeric(df[x], errors="coerce")
        y_num = pd.to_numeric(df[y], errors="coerce")
        m = x_num.notna() & y_num.notna()
        if m.sum() > 2:
            b1, b0 = np.polyfit(x_num[m], y_num[m], 1)
            xs = np.linspace(x_num[m].min(), x_num[m].max(), 50)
            fig.add_trace(go.Scatter(x=xs, y=b1 * xs + b0, mode="lines",
                                     name="tendencia", line=dict(color=PALETTE["coral"], dash="dot")))
    return _apply(fig, title, ytitle, xtitle, height)


def line(df, x, y, title=None, color=None, ytitle="", xtitle="", height=440, markers=True):
    fig = px.line(df, x=x, y=y, color=color, markers=markers,
                  color_discrete_sequence=CHART_SEQUENCE)
    return _apply(fig, title, ytitle, xtitle, height)


def heatmap(matrix, title=None, zmin=None, zmax=None, labels=None, height=460,
            colorscale="Teal"):
    fig = px.imshow(matrix, text_auto=".2f", color_continuous_scale=colorscale,
                    zmin=zmin, zmax=zmax, aspect="auto")
    if labels:
        fig.update_xaxes(title_text=labels[0])
        fig.update_yaxes(title_text=labels[1])
    return _apply(fig, title, height=height)


def radar(categories, values, title=None, height=440, color=None, fill=True):
    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(r=values + values[:1], theta=categories + categories[:1],
                                  fill="toself" if fill else None,
                                  line=dict(color=color or PALETTE["teal"], width=2.5),
                                  name="perfil"))
    fig.update_layout(polar=dict(bgcolor="white",
                                 radialaxis=dict(gridcolor=PALETTE["line"], range=[0, 100] if max(values, default=100) > 1.5 else [0, max(values) * 1.2])))
    return _apply(fig, title, height=height)


def gauge(value, title=None, delta_label=None, height=260, thresholds=None, value_range=(0, 100)):
    thresholds = thresholds or [(0, "#FF6B6B"), (30, "#F08A24"), (45, "#E9B44C"),
                                (60, "#00A896"), (80, "#3E9B4F")]
    fig = go.Figure(go.Indicator(
        mode="gauge+number", value=value,
        title=dict(text=title or "", font=dict(size=14)),
        number=dict(valueformat=".1f"),
        gauge=dict(axis=dict(range=value_range),
                   bar=dict(color=PALETTE["navy"], thickness=.28),
                   steps=[dict(range=[thresholds[i][0], thresholds[i + 1][0]],
                               color=thresholds[i][1])
                          for i in range(len(thresholds) - 1)])))
    fig.update_layout(height=height, margin=dict(l=24, r=24, t=70, b=24))
    return _apply(fig)


def stacked(df, x, y, color, title=None, ytitle="", xtitle="", height=440, normalize=False):
    d = df.copy()
    if normalize:
        d[y] = d[y] / d.groupby(x)[y].transform("sum") * 100
    fig = px.bar(d, x=x, y=y, color=color, barmode="stack",
                 color_discrete_sequence=CHART_SEQUENCE)
    return _apply(fig, title, ytitle, xtitle, height)


def crosstab_heat(df, index, columns, values, aggfunc="count", title=None, height=460):
    ct = pd.crosstab(df[index], df[columns], values=df[values], aggfunc=aggfunc)
    if ct.empty:
        ct = pd.crosstab(df[index], df[columns])
    return heatmap(ct, title=title, height=height)
