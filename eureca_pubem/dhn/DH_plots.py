import matplotlib.pyplot as plt
import numpy as np


def plot_dhn_capacity_layout(
    subsystems: dict,
    *,
    streets_gdf=None,
    unassigned_nodes=None,
    figsize=(12, 12),
    save_path=None,
):
    """
    Plot DHN subsystems with full visual enhancements.

    Parameters
    ----------
    subsystems : dict[str, DistrictHeating]
        One subsystem per supply.
    streets_gdf : GeoDataFrame, optional
        Street geometry for background.
    unassigned_nodes : list[Node], optional
        Consumers not assigned to any supply.
    figsize : tuple
        Figure size.
    save_path : str, optional
        If given, saves figure to file.
    """

    fig, ax = plt.subplots(figsize=figsize)

    # background streets
    if streets_gdf is not None:
        streets_gdf.plot(
            ax=ax,
            color="lightgrey",
            linewidth=0.6,
            alpha=0.35,
            zorder=0,
        )

    cmap = plt.cm.get_cmap("tab10")

    legend_used = set()

    for i, (name, system) in enumerate(subsystems.items()):
        color = cmap(i % 10)

        # ---- plot lines ----
        for l in system.lines:
            n1 = next(n for n in system.nodes if n.node_id == l.start_node)
            n2 = next(n for n in system.nodes if n.node_id == l.end_node)

            ax.plot(
                [n1.x, n2.x],
                [n1.y, n2.y],
                color=color,
                linewidth=1.5,
                alpha=0.85,
                zorder=1,
            )

        # ---- plot nodes ----
        for n in system.nodes:
            if n.node_type == "supply":
                size = 120
                if hasattr(n, "capacity"):
                    size = 120 + 0.03 * n.capacity

                label = name if name not in legend_used else None
                legend_used.add(name)

                ax.scatter(
                    n.x,
                    n.y,
                    s=size,
                    color=color,
                    edgecolor="black",
                    linewidth=1.2,
                    zorder=5,
                    label=label,
                )

            elif n.node_type == "consumer":
                size = 25
                if hasattr(n, "peak_demand"):
                    size = 15 + 0.06 * n.peak_demand

                ax.scatter(
                    n.x,
                    n.y,
                    s=size,
                    color=color,
                    alpha=0.85,
                    zorder=3,
                )

    # ---- unassigned consumers ----
    if unassigned_nodes:
        ax.scatter(
            [n.x for n in unassigned_nodes],
            [n.y for n in unassigned_nodes],
            s=30,
            color="grey",
            alpha=0.6,
            zorder=2,
            label="Unassigned consumers",
        )

    # ---- aesthetics ----
    ax.set_aspect("equal")
    ax.set_title("District Heating Network – Capacity-Based Assignment", fontsize=14)
    ax.axis("off")
    ax.legend(frameon=True, loc="best")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)

    plt.show()
