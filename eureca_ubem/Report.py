from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.units import cm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import Paragraph, Frame, Spacer, Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.utils import ImageReader
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
from io import BytesIO


def create_scatter_image(x, y, x_label, y_label, title):
    """
    Create an in-memory scatter plot and return it as an ImageReader object.
    """
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    ax.scatter(x, y, color='teal', alpha=0.8)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(True)
    buf = BytesIO()
    canvas = FigureCanvas(fig)
    canvas.print_png(buf)
    plt.close(fig)
    buf.seek(0)
    return ImageReader(buf)


def generate_summary_pdf(output_path, scenario_dict, summary_dict):
    """
    Generate a PDF report summarizing each scenario with descriptive text,
    a statistical summary table, and 3 inline scatter plots.

    Parameters
    ----------
    output_path : str
        File path for saving the PDF.
    scenario_dict : dict
        Original intervention definitions per scenario.
    summary_dict : dict
        Calculated metrics per scenario (from summarize_scenario_outputs).
    """
    c = canvas.Canvas(output_path, pagesize=A4)
    width, height = A4
    styles = getSampleStyleSheet()
    normal_style = styles["Normal"]
    bold_style = ParagraphStyle(name="Bold", parent=normal_style, fontName="Helvetica-Bold")

    margin = 2 * cm
    frame_width = width - 2 * margin

    for scenario_name in scenario_dict:
        interventions = scenario_dict[scenario_name]
        summary = summary_dict[scenario_name]

        # TITLE
        c.setFont("Helvetica-Bold", 18)
        c.drawString(margin, height - margin, f"Scenario: {scenario_name.replace('_', ' ').title()}")

        # Build description text
        description = "In this scenario:"
        for k, v in interventions.items():
            if k.lower() == "pv":
                description += f"<br/>• About {v}% of buildings have photovoltaics installed on their rooftops."
            elif k.lower() == "ts":
                description += f"<br/>• Around {v}% of buildings have thermal solar panels."
            elif k.lower() == "hvac":
                description += f"<br/>• {v}% of buildings received HVAC system upgrades."
            elif k.lower() == "envelope":
                description += f"<br/>• {v}% of buildings underwent envelope insulation improvements."
            elif k.lower() == "deep_retrofit":
                description += f"<br/>• {v}% of buildings received deep retrofit interventions."
            else:
                description += f"<br/>• {v}% of buildings include the intervention: {k}."

        # Description Frame
        description_paragraph = Paragraph(description, style=normal_style)
        description_frame = Frame(margin, height - 7*cm, frame_width, 5*cm, showBoundary=0)
        description_frame.addFromList([description_paragraph], c)

        # Metrics Title
        c.setFont("Helvetica-Bold", 13)
        c.drawString(margin, height - 7.5*cm, "Key Performance Metrics:")

        # Statistical summary table
        metrics_data = [
            ["Metric", "Mean", "Max", "Min", "Std. Dev"],
            ["Cost (€)",
             f"{summary['cost_mean']:.2f}",
             f"{summary['cost_max']:.2f}",
             f"{summary['cost_min']:.2f}",
             f"{summary['cost_std']:.2f}"],
            ["Primary Energy (kWh)",
             f"{summary['primary_energy_mean']:.2f}",
             f"{summary['primary_energy_max']:.2f}",
             f"{summary['primary_energy_min']:.2f}",
             f"{summary['primary_energy_std']:.2f}"],
            ["Non-Renewable PE (kWh)",
             f"{summary['primary_energy_nren_mean']:.2f}",
             f"{summary['primary_energy_nren_max']:.2f}",
             f"{summary['primary_energy_nren_min']:.2f}",
             f"{summary['primary_energy_nren_std']:.2f}"],
            ["CO₂ Emission (tons)",
             f"{summary['co2_emission_mean']:.2f}",
             f"{summary['co2_emission_max']:.2f}",
             f"{summary['co2_emission_min']:.2f}",
             f"{summary['co2_emission_std']:.2f}"],
        ]

        table = Table(metrics_data, colWidths=[6*cm, 3*cm, 3*cm, 3*cm, 3*cm])
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.grey),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
            ("ALIGN", (1, 1), (-1, -1), "CENTER"),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 6),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
        ]))

        # Position metrics table just after the description
        metrics_frame = Frame(margin, height - 15.5*cm, frame_width, 6.5*cm, showBoundary=0)
        metrics_frame.addFromList([Spacer(1, 0.3*cm), table], c)

        # Inline scatter plots
        cost = summary['cost']
        scatter_y_list = [
            ("CO₂ Emissions (tons)", summary['co2_emission'], "CO₂ vs Cost"),
            ("Primary Energy (kWh)", summary['primary_energy'], "Primary Energy vs Cost"),
            ("Non-Renewable PE (kWh)", summary['primary_energy_nren'], "PE_nren vs Cost")
        ]

        chart_y = 6.2 * cm
        chart_width = (frame_width - 2 * cm) / 3
        chart_height = 3.5 * cm

        for i, (y_label, y_data, title) in enumerate(scatter_y_list):
            img = create_scatter_image(cost, y_data, "Cost (€)", y_label, title)
            c.drawImage(img,
                        margin + i * (chart_width + 1 * cm),
                        chart_y,
                        width=chart_width,
                        height=chart_height,
                        preserveAspectRatio=True)

        c.showPage()

    c.save()
    return output_path
