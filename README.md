# Clustering Code for Pixel Detectors

## Description

### Overview
This repository contains clustering code for pixel detectors used in particle physics. The primary goal is to implement algorithms for clustering pixels in detector data, improving the analysis and understanding of particle interactions and energy deposition. This code can be used for applications in radiotherapy, where precise knowledge of energy spectra is crucial.

### Features and Applications
- Efficient clustering of pixel data.
- Applicable in particle physics and radiotherapy.
- Supports multiple simultaneous frame clustering.
- Energy calibration implementation.
- Output in both ROOT files and tabular text format.
- Generation of histograms and other descriptive graphical analysis.

For a detailed description of the methodology and application, refer to [link](https://www.if.ufrj.br/wp-content/uploads/2020/11/Dissertacao_Marcos_Vieira_IF_FisAplic_UFRJ_final_ficha.pdf).

## Repository Structure

### Key Files

- **clustering.cpp**: Main clustering function without calibration or multi-frame support.
- **clustering_cal.cpp**: Clustering function with calibration support.
- **clustering_write_txt.cpp**: Clustering function with text output.
- **clustering_multiframes.cpp**: Clustering function supporting multiple frames.
- **clusteringmacro.C**: Macro file to run the clustering code.

### Directories
- **input/**: Directory containing input data files, e.g. a `.tpx` file.
- **output/**: Directory where the output files, such as ROOT and text files, are saved.

## Usage

### Input
The input data consists of pixel detector readings, e.g. a `.tpx` file.

### Output
The output is a set of clustered pixel data, which can be further analyzed or visualized. The output format varies depending on the selected clustering function, with options for ROOT and text files.

### Example
To use the clustering code, follow these steps:

1. **Clone the repository**

2. **Navigate to the repository**

3. **Ensure you have [ROOT](https://root.cern/) installed and configured.**

4. **Unzip the example `input/cs137.zip`:** This will extract `cs137.tpx` file to use in the code.

5. **Comment/Uncomment the selected code in the macro `clusteringmacro.C`** to choose the desired clustering code.

6. **Run the `clusteringmacro.C` macro:**
    ```sh
    root -l -b -q clusteringmacro.C
    ```

### Requirements
- **ROOT Framework**: [ROOT Installation Guide](https://root.cern/install/)

  - The code has been tested under [ROOT Release 6.32/00 - 2024-05-28](https://root.cern/releases/release-63200/).
  - The dependencies are presented [here.](https://root.cern/install/dependencies/)

- **Operating System**: The code has been tested on Windows 10 x64 and Linux distributions.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contribution
Contributions are welcome! Please follow these steps:
1. Fork the repository.
2. Create a new branch:
    ```sh
    git checkout -b feature-branch
    ```
3. Commit your changes:
    ```sh
    git commit -am 'Add new feature'
    ```
4. Push to the branch:
    ```sh
    git push origin feature-branch
    ```
5. Create a new Pull Request.

## Authors
- [Marcos Vieira](https://github.com/marcosgvjunior)

## References
- [ROOT Data Analysis Framework](https://root.cern)
- [Marcos Vieira's Dissertation](https://www.if.ufrj.br/wp-content/uploads/2020/11/Dissertacao_Marcos_Vieira_IF_FisAplic_UFRJ_final_ficha.pdf)
