import pandas as pd
import matplotlib.pyplot as plt

def blank_plot(message):
    """
    Create a blank plot with a centered text message.

    Parameters:
        message (str): The message to be displayed in the plot.

    Returns:
        matplotlib.axes._subplots.AxesSubplot: The blank plot with the message.
    """

    # Create a DataFrame with a single row containing the message
    df = pd.DataFrame({"a": [0], "b": [0], "n": [message]})

    # Create a new figure
    plt.figure()

    # Add the text to the plot at position (0, 0)
    plt.text(0, 0, message, ha="center", va="center")

    # Remove the labels from x-axis and y-axis
    plt.xlabel(None)
    plt.ylabel(None)

    # Remove the ticks from x-axis and y-axis
    plt.xticks([])
    plt.yticks([])

    # Adjust the plot layout
    plt.tight_layout()

    # Return the current axes
    return plt.gca()
