from napari.utils.notifications import show_error
from magicgui import magic_factory
import matplotlib.pyplot as plt
import numpy as np
import napari

def bbox(rectangle):
    if rectangle.shape != (4, 2):
        show_error("Please, define rectangular shapes!")
        return None
    else:
        # Extract cropping indices from the shape information
        row_min = int(rectangle[:, 0].round().min())
        row_max = int(rectangle[:, 0].round().max())
        col_min = int(rectangle[:, 1].round().min())
        col_max = int(rectangle[:, 1].round().max())
        return (row_min, row_max, col_min, col_max)
    
def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text

@magic_factory(
    call_button='Calculate',
    Display_Mode={"choices": ['Amplitude', 'Phase', 'Amplitude + Phase']},
)
def fft2(
    Amplitude_Image: 'napari.layers.image.Image',
    Phase_Image: 'napari.layers.image.Image',
    Selected_Zones: 'napari.layers.shapes.Shapes',
    Display_Mode: str = 'Amplitude + Phase'
) -> None:
    """
    Compute the 2D FFT of the selected region in an image.

    Parameters
    ----------
    Amplitude_Image : napari.layers.image.Image
        Magnitude information of the image.
    Phase_Image : napari.layers.image.Image
        Phase information of the image.
    Selected_Zones : napari.layers.shapes.Shapes
        The selected regions to compute the 2D FFT.
    Display_Mode : str, optional
        Display mode for the FFT ('Amplitude', 'Phase', or 'Amplitude + Phase').

    Returns
    -------
    None
    """
    for shape in Selected_Zones.data:
        # Extract cropping indices
        indices = bbox(shape)

        if indices:
            # Validate input shapes and indices
            if not np.all(Amplitude_Image.data[0].shape == Phase_Image.data[0].shape):
                show_error('The amplitude image and the phase image must have the same size')
                return
            elif (not 0 <= indices[0] < indices[1] <= Amplitude_Image.data[0].shape[0]) or (not 0 <= indices[2] < indices[3] <= Amplitude_Image.data[0].shape[1]):
                show_error('The selected region must be included in the image')
                return

            # Crop the data and apply phase information
            data_cropped = Amplitude_Image.data[0][indices[0]:indices[1], indices[2]:indices[3]] * np.exp(1j * Phase_Image.data[0][indices[0]:indices[1], indices[2]:indices[3]])

            # Perform 2D FFT and center the spectrum if required
            data_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(data_cropped)))

            # Display the result based on the selected mode
            amplitude = np.abs(data_fft)
            phase = np.angle(data_fft)
            if Display_Mode == 'Amplitude':
                data_fft_display = amplitude

            elif Display_Mode == 'Phase':
                data_fft_display = phase

            elif Display_Mode == 'Amplitude + Phase': # HSV representation
                phase_normalized = (phase + np.pi) / (2 * np.pi)
                amplitude_normalized = amplitude / min(np.max(amplitude), np.mean(amplitude) + 3 * np.std(amplitude))
                np.clip(amplitude_normalized, 0, 1, out=amplitude_normalized)
                data_fft_display = plt.matplotlib.colors.hsv_to_rgb(
                    np.stack([phase_normalized, np.ones_like(amplitude_normalized), amplitude_normalized], axis=-1)
                )

            else:
                raise ValueError("Invalid Display_Mode. Use 'Amplitude', 'Phase', or 'Amplitude + Phase'.")

            # Display the result using matplotlib
            filename = remove_prefix(Amplitude_Image.name, '[abs] ')
            _, ax = plt.subplots()
            ax.set_title(f"{filename}\nFTT2D ({Display_Mode})")
            ax.imshow(data_fft_display, extent=[-0.5, 0.5, -0.5, 0.5])
            ax.set_xlabel("Normalized Frequency")
            ax.set_ylabel("Normalized Frequency")
            ax.axis("scaled")
            ax.axis([-0.5, 0.5, -0.5, 0.5])
            plt.show()

@magic_factory(
    call_button='Plot',
)
def hist(
    Amplitude_Image: 'napari.layers.image.Image',
    Selected_Zones: 'napari.layers.shapes.Shapes',
    Bins: int = 1024
) -> None:
    """
    Plot histogram of the selected region in an image.

    Parameters
    ----------
    Amplitude_Image : napari.layers.image.Image
        Magnitude information of the image.
    Selected_Zones : napari.layers.shapes.Shapes
        The selected regions to compute the histogram.
    Bins : int, optional
        Number of Bins in the histogram, by default 1024.

    Returns
    -------
    None
        The function displays the histogram plot.
    """
    for shape in Selected_Zones.data:
        # Extract cropping indices
        indices = bbox(shape)

        if indices:
            # Validate indices
            if (not 0 <= indices[0] < indices[1] <= Amplitude_Image.data[0].shape[0]) or (not 0 <= indices[2] < indices[3] <= Amplitude_Image.data[0].shape[1]):
                show_error('The selected region must be included in the image')
                return

            # Crop the data
            data_cropped = Amplitude_Image.data[0][indices[0]:indices[1], indices[2]:indices[3]]

            # Display the result using matplotlib
            _, ax = plt.subplots()
            ax.set_title('Histogram')
            ax.hist(
                data_cropped.flatten(),
                bins=Bins,
                range=(data_cropped.min(), data_cropped.max()),
                histtype="stepfilled",
                edgecolor="none",
                density=True
            )
            ax.set_xlabel('Values')
            ax.set_ylabel('Number of pixels')
            plt.show()