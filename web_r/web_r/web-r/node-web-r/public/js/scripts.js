import {WebR} from 'https://webr.r-wasm.org/latest/webr.mjs';
import * as Plot from "https://cdn.jsdelivr.net/npm/@observablehq/plot@0.6/+esm";

const webR = new WebR();
await webR.init();

let exploreMotifScript = () => {
    return '\n' +
        'filter.data.frame <- function(data, col.number.to.filter, filter.value){\n' +
        '  filtered_df <- subset(data, data[, col.number.to.filter] == filter.value)\n' +
        '  return(filtered_df)\n' +
        '}\n' +
        '\n' +
        '\n' +
        'time.series.from.dataframe <- function(data, col.number, col.number.end, row.start, row.end){\n' +
        '  if(col.number.end != -1){\n' +
        '    time.series.from.dataframe.multu.column(data, col.number, col.number.end, row.start, row.end)\n' +
        '  }else {\n' +
        '    df <- data.frame(matrix(ncol = 1, nrow = 0))\n' +
        '    rows <- nrow(data)\n' +
        '    \n' +
        '    for (i in 1:rows) {\n' +
        '      UNRATE <- data[i, col.number]\n' +
        '      df <- rbind(df, data.frame(UNRATE=UNRATE))\n' +
        '    }\n' +
        '    \n' +
        '    return(df)\n' +
        '  }\n' +
        '}\n' +
        '\n' +
        'time.series.from.dataframe.multu.column <- function(data, col.number, col.number.end, row.start, row.end){\n' +
        '  df <- data.frame(matrix(ncol = 1, nrow = 0))  \n' +
        '  # rows <- nrow(data)\n' +
        '  # \n' +
        '  # for (i in 1:rows) {\n' +
        '  #   # chunk <- data[1:12, i]\n' +
        '  #   UNRATE <- data[i, col.number]\n' +
        '  #   df <- rbind(df, data.frame(UNRATE=UNRATE))\n' +
        '  # }\n' +
        '  \n' +
        '  # columns <- ncol(data)\n' +
        '  \n' +
        '  for (i in col.number:col.number.end) {  # 2:columns\n' +
        '    chunk <- data[row.start:row.end, i]  # 1:12\n' +
        '    for (j in row.start:row.end){  # 1:12\n' +
        '      # print(chunk[j])\n' +
        '      UNRATE <- chunk[j]\n' +
        '      df <- rbind(df, data.frame(UNRATE=UNRATE))\n' +
        '    }\n' +
        '  }\n' +
        '  \n' +
        '  return(df)\n' +
        '}\n' +
        '\n' +
        '\n' +
        'euclidean <- function(a, b) sqrt(sum((a - b)^2))\n' +
        '\n' +
        '\n' +
        'correlation.distance <- function(p, q) {\n' +
        '  diff.p <- diff(p)\n' +
        '  diff.q <- diff(q)\n' +
        '  pq.production <- diff.p * diff.q\n' +
        '  p.square <- diff.p * diff.p\n' +
        '  q.square <- diff.q * diff.q \n' +
        '  s1 <- 0\n' +
        '  s2 <- 0\n' +
        '  s3 <- 0\n' +
        '  n <- length(p)\n' +
        '  for(i in 1: (n - 1)){\n' +
        '    s1 <- s1 + pq.production[i]\n' +
        '    s2 <- s2 + p.square[i]\n' +
        '    s3 <- s3 + q.square[i]\n' +
        '  }\n' +
        '  \n' +
        '  s <- s2 * s3\n' +
        '  s <- sqrt(s)\n' +
        '  \n' +
        '  return(s1 / s)\n' +
        '}\n' +
        '\n' +
        '\n' +
        'chouakria.similarity.measure <- function(p, q, k) {\n' +
        '  cd <- correlation.distance(p, q)\n' +
        '  ch <- 2 / (1 + exp((k * cd)))\n' +
        '  euc <- euclidean(p, q)\n' +
        '  \n' +
        '  euc * ch\n' +
        '  return(euc * ch)\n' +
        '}\n' +
        '\n' +
        '\n' +
        'distance <- function(p, q, dstance.name) {\n' +
        '  if(dstance.name == "euclidean"){\n' +
        '    return(euclidean(p, q))\n' +
        '  }\n' +
        '  else if(dstance.name == "chouakria"){\n' +
        '    k <- 1\n' +
        '    return(chouakria.similarity.measure(p, q, k))\n' +
        '  } else {\n' +
        '    return(-1)\n' +
        '  }\n' +
        '}\n' +
        '\n' +
        '\n' +
        'generateSimilarityMatrix <- function(sequencLength, df, exclusion_zone, dstance.name) {\n' +
        '  n <-  (nrow(df) - sequencLength) + 1\n' +
        '  distance.matrix <- matrix(0, n, n)\n' +
        '  \n' +
        '  for (j in 1 : (nrow(df) - sequencLength)) {\n' +
        '    timeSeries.sequence <- df[(j) : (j + (sequencLength - 1)) , 1]\n' +
        '    exclusion.interval <- j + exclusion_zone\n' +
        '    \n' +
        '    for (i in j : ((nrow(df) - sequencLength) + 1)){\n' +
        '      chunk <- df[i : (i + (sequencLength - 1)) , 1]\n' +
        '      if(i > exclusion.interval){\n' +
        '        # distance <- euclidean(chunk, timeSeries.sequence)\n' +
        '        distance <- distance(chunk, timeSeries.sequence, dstance.name)\n' +
        '        distance.matrix[j, i] <- distance\n' +
        '        # distance.matrix[i, j] <- distance\n' +
        '      }\n' +
        '    }\n' +
        '  }\n' +
        '  \n' +
        '  return (distance.matrix)\n' +
        '}\n' +
        '\n' +
        '\n' +
        'miniGreaterThanZero <- function(similarity.matrix, row.index, threshold) {\n' +
        '  n <- ncol(similarity.matrix)\n' +
        '  distance.vector <- replicate(n, 0)\n' +
        '  \n' +
        '  non_zero_elements <- similarity.matrix[row.index, ( similarity.matrix[row.index, ] != 0 & similarity.matrix[row.index, ] <= threshold ) ]\n' +
        '  \n' +
        '  if(length(non_zero_elements) == 0 ){\n' +
        '    return(distance.vector)    \n' +
        '  }\n' +
        '  \n' +
        '  for(j in 1 : length(non_zero_elements) ){\n' +
        '    for(i in row.index : n ){\n' +
        '      if( !is.na(similarity.matrix[row.index, i]) & !is.na(non_zero_elements[j]) & similarity.matrix[row.index, i] == non_zero_elements[j] ){\n' +
        '        distance.vector <- replace(distance.vector, i, i)\n' +
        '      }\n' +
        '    }\n' +
        '  }    \n' +
        '  return(distance.vector)\n' +
        '}\n' +
        '\n' +
        'generate.plots <- function(similarity.matrix, epsilon.query, k, n, ylabel, mainTitle, df, sequence.length) {\n' +
        '  \n' +
        ' print(nrow(similarity.matrix))\n' +
        '  n <- min(n, nrow(similarity.matrix))\n' +
        '  for(i in k : n){\n' +
        '    distance.vector <- miniGreaterThanZero(similarity.matrix, i, epsilon.query)\n' +
        '    epsilon <- distance.vector[distance.vector != 0]\n' +
        '    \n' +
        '    if(length(epsilon) > 0){\n' +
        '      similarity.number <- length(epsilon)\n' +
        '      \n' +
        '      par(mfrow = c(2, 2))\n' +
        '      layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE))\n' +
        '      \n' +
        '      # var = readline(prompt = "press yes to acept and false to refuse : ");\n' +
        '      # print(var)\n' +
        '      \n' +
        '      # plot(c(1:length(df[,1])), df[,1], type = "l", xlab = "Index", ylab = "Values", col = "black", main = "Multiple Time Series")\n' +
        '      plot(c(1:length(df[,1])), df[,1], type = "l", xlab = "Time Index", ylab = ylabel, col = "black", main = mainTitle)\n' +
        '      lines(c(i : (i + sequence.length - 1)), df[i : (i + sequence.length - 1), 1], col = "green")\n' +
        '      for(j in 1: similarity.number){\n' +
        '        lines(c(epsilon[j] : (epsilon[j] + sequence.length - 1)), df[epsilon[j] : (epsilon[j] + sequence.length - 1), 1], col = "red")\n' +
        '      }\n' +
        '      \n' +
        '      plot(as.ts(df[i : (i + sequence.length - 1), 1]), type = "l", xlab = "Date", ylab = "Values", col = "black", main = "Motifs")\n' +
        '      for(j in 1: similarity.number){\n' +
        '        lines(as.ts(df[epsilon[j] : (epsilon[j] + sequence.length - 1), 1]), col = "blue")\n' +
        '      }\n' +
        '    } \n' +
        '  }\n' +
        '}\n' +
        '\n' +
        '\n' +
        'normalize.timeseries <- function(df) {\n' +
        '  m <- mean(df[ , 1])\n' +
        '  std <- sd(df[ , 1])\n' +
        '  df[ , 1] <- df[ , 1] - m\n' +
        '  df[ , 1] <- df[ , 1] / std\n' +
        '  \n' +
        '  return(df)\n' +
        '}\n' +
        '\n' +
        '\n' +
        'explore.motif <- function(data.path, column.number.to.filter, filter.value, column.number.to.extract.data, sequence.length, distance.name, epsilon.query, k, n, col.number.end, row.start, row.end) {\n' +
        '  df <- get.timeseries(data.path, column.number.to.filter, filter.value, column.number.to.extract.data, col.number.end, row.start, row.end)\n' +
        '  similarity.matrix <- generateSimilarityMatrix(sequence.length, df, sequence.length, distance.name)  # euclidean # chouakria\n' +
        '  \n' +
        '  # return(similarity.matrix)\n' +
        '  my_dataframe <- as.data.frame(similarity.matrix)\n' +
        '  return(my_dataframe)\n' +
        '}\n' +
        '\n' +
        '# get.timeseries <- function(data.path, column.number.to.filter, filter.value, column.number.to.extract.data) {\n' +
        '#   data <- read.csv(data.path)\n' +
        '#   \n' +
        '#   filtered_df <- filter.data.frame(data, column.number.to.filter, filter.value)\n' +
        '#   df <- time.series.from.dataframe(filtered_df, column.number.to.extract.data)\n' +
        '#   df <- normalize.timeseries(df)\n' +
        '#   \n' +
        '#   return(df)\n' +
        '# }\n' +
        '\n' +
        'get.timeseries <- function(data.path, column.number.to.filter, filter.value, column.number.to.extract.data, col.number.end, row.start, row.end) {\n' +
        '  data <- read.csv(data.path)\n' +
        '  \n' +
        '  if(column.number.to.filter != -1){\n' +
        '    filtered_df <- filter.data.frame(data, column.number.to.filter, filter.value)\n' +
        '    # df <- time.series.from.dataframe(filtered_df, column.number.to.extract.data) \n' +
        '    df <- time.series.from.dataframe(filtered_df, column.number.to.extract.data, col.number.end, row.start, row.end)\n' +
        '  }else{\n' +
        '    # df <- time.series.from.dataframe(data, column.number.to.extract.data)\n' +
        '    df <- time.series.from.dataframe(data, column.number.to.extract.data, col.number.end, row.start, row.end)\n' +
        '  }\n' +
        '  \n' +
        '  df <- normalize.timeseries(df)\n' +
        '  \n' +
        '  return(df)\n' +
        '}\n' +
        '\n' +
        'get.timeseries.data <- function(data.path, column.number.to.filter, filter.value, column.number.to.extract.data, col.number.end, row.start, row.end) {\n' +
        '  data <- get.timeseries(data.path, column.number.to.filter, filter.value, column.number.to.extract.data, col.number.end, row.start, row.end)\n' +
        '    return(data[,1]) \n' +
        '}\n' +
        '\n' +
        'get.data.bruto <- function(data.path){\n' +
        '  data <- read.csv(data.path)\n' +
        '  return(data)\n' +
        '}\n';
}

// let rscript = await exploreMotifScript();
// await webR.evalR(rscript);
//
// let returnValue = await webR.evalR('explore.motif("http://localhost:3000/csv/annual-co2-emissions-per-country", 2, "ALB", 4, 10, "euclidean", 0.5, 1, 100)')
// let x = await returnValue.toArray();

function arrayToMatrix(arr, rows, cols) {
    const matrix = [];
    let k = 0;
    for (let i = 0; i < rows; i++) {
        matrix.push([]);
    }
    for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
            matrix[j][i] = arr[k];
            k++;
        }
    }
    return matrix;
}

const grafic = (data, selectorString, i, distanceArray, sequenceLength) => {
    let graphMetadata = {
        title: 'Lindjet 1990 - 2022',
        subTitle: 'Lindjet per cdo muaj',
        xlable: 'indexes',
        ylable: 'Lindje'
    }


    // Set the dimensions and margins of the graph
    const margin = {top: 70, right: 50, bottom: 50, left: 50},
        width = 1200 - margin.left - margin.right,
        height = 500 - margin.top - margin.bottom;

    // Append the SVG object to the body of the page
    const svg = d3.select(selectorString)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

    let mainSequenceDatas = data.filter(dataElement => dataElement.index >= i && dataElement.index < (i + sequenceLength))
    let patternDataArrays = [];

    for (let distanceIndex of distanceArray) {
        let dataArray = data.filter(dataElement => dataElement.index >= distanceIndex && dataElement.index < (distanceIndex + sequenceLength))
        patternDataArrays.push(dataArray);
    }

    // X scale
    const x = d3.scaleTime()
        // .domain(d3.extent(data, d => d.date))
        .domain(d3.extent(data, d => d.index))
        .range([0, width]);

    // Y scale
    const y = d3.scaleLinear()
        // .domain([0, d3.max(data, d => d.value)])
        .domain([d3.min(data, d => d.value), d3.max(data, d => d.value)])
        .range([height, 0]);

    // Define the line
    const line = d3.line()
        // .x(d => x(d.date))
        .x(d => x(d.index))
        .y(d => y(d.value));

    // Add the X Axis
    svg.append("g")
        .attr("transform", `translate(0,${height})`)
        .call(d3.axisBottom(x)
            // .ticks((data.length * 0.5)).tickFormat(d3.format("d")));
            .ticks((data.length * 0.2)).tickFormat(d3.format("d")));

    // Add the Y Axis
    svg.append("g")
        .call(d3.axisLeft(y));

    // Add the line path
    svg.append("path")
        .datum(data)
        .attr("class", "line")
        .attr("d", line)
        .style("stroke", "black");

    // Add the line path
    svg.append("path")
        .datum(mainSequenceDatas)
        .attr("class", "line")
        .attr("d", line)
        .style("stroke", "red");

    for (let arr of patternDataArrays) {
        // Add the line path
        svg.append("path")
            .datum(arr)
            .attr("class", "line")
            .attr("d", line)
            .style("stroke", "green");
    }

    // Add chart title
    svg.append("text")
        .attr("x", width / 2)
        .attr("y", ((-margin.top / 2) - 15))
        .attr("class", "title")
        .text(graphMetadata.title);

    // Add chart subtitle
    svg.append("text")
        .attr("x", width / 2)
        .attr("y", -margin.top / 2 + 10)
        .attr("class", "subtitle")
        .text(graphMetadata.subTitle);

    // Add legend
    const legend = svg.append("g")
        .attr("class", "legend")
        .attr("transform", `translate(${width - 200},${margin.top - margin.top})`);

    legend.append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", 10)
        .attr("height", 10)
        .style("fill", "green");

    legend.append("text")
        .attr("x", 20)
        .attr("y", 10)
        .text("Sequence that match motifs");

    legend.append("rect")
        .attr("x", 0)
        .attr("y", 20)
        .attr("width", 10)
        .attr("height", 10)
        .style("fill", "red");

    legend.append("text")
        .attr("x", 20)
        .attr("y", 30)
        .text("Found motifs");

    // x label
    svg.append("text")
        .attr("x", width / 2)
        .attr("y", height + 50)
        .attr("text-anchor", "middle")
        .attr("transform", `rotate(0, ${width / 2}, ${height + 30})`)
        .text(graphMetadata.xlable);

    // y label
    svg.append("text")
        .attr("x", -30)
        .attr("y", (height / 2) - 30)
        .attr("text-anchor", "middle")
        .attr("transform", `rotate(-90, ${-10}, ${height / 2})`)
        .text(graphMetadata.ylable);

}

const patternVizualzer = (data, selectorString, i, distanceArray, sequenceLength) => {

    // Set the dimensions and margins of the graph
    const margin = {top: 40, right: 30, bottom: 40, left: 50},
        width = 300 - margin.left - margin.right,
        height = 300 - margin.top - margin.bottom;

    // Append the SVG object to the body of the page
    const svg = d3.select(selectorString)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

    let patternDataArray = [];
    let mainSequenceData = data.filter(dataElement => dataElement.index >= i && dataElement.index < (i + sequenceLength))
    mainSequenceData = mainSequenceData.map((dataElement, index) => {
        let temp = Object.assign({}, dataElement);
        temp.index = index;
        return temp;
    });

    // Defining colors for charts
    const customColors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
        "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
    ];
    const colorScale = d3.scaleOrdinal(customColors);

    let j = 1;
    for (let distanceIndex of distanceArray) {
        let dataArray = data.filter(dataElement => dataElement.index >= distanceIndex && dataElement.index < (distanceIndex + sequenceLength))
        dataArray = dataArray.map((dataElement, index) => {
            let temp = Object.assign({}, dataElement);
            temp.index = index;
            return temp;
        });
        patternDataArray.push(dataArray);
    }

    // X scale
    const x = d3.scaleTime()
        .domain(d3.extent(mainSequenceData, d => d.index))
        .range([0, width]);

    // Y scale
    const y = d3.scaleLinear()
        .domain([d3.min(mainSequenceData, d => d.value), d3.max(mainSequenceData, d => d.value)])
        .range([height, 0]);

    // Define the line
    const line = d3.line()
        .x(d => x(d.index))
        .y(d => y(d.value));

    // Add the X Axis
    svg.append("g")
        .attr("transform", `translate(0,${height})`)
        .call(d3.axisBottom(x)
            .ticks((mainSequenceData.length * 0.3)).tickFormat(d3.format("d")));

    // Add the Y Axis
    svg.append("g")
        .call(d3.axisLeft(y));

    // Add the line path
    svg.append("path")
        .datum(mainSequenceData)
        .attr("class", "line")
        .attr("d", line)
        .style("stroke", "black");

    for (let arr of patternDataArray) {
        let xMainSequence = d3.scaleTime()
            .domain(d3.extent(arr, d => d.index))
            .range([0, width]);

        // Y scale
        let yMainSequence = d3.scaleLinear()
            .domain([d3.min(arr, d => d.value), d3.max(arr, d => d.value)])
            .range([height, 0]);

        // Define the line
        let mainSequenceLine = d3.line()
            .x(d => xMainSequence(d.index))
            .y(d => yMainSequence(d.value));

        // Add the line path
        svg.append("path")
            .datum(arr)
            .attr("class", "line")
            .attr("d", mainSequenceLine)
            .style("stroke", colorScale(j));
        j++;
    }

    // Add chart title
    svg.append("text")
        .attr("x", width / 2)
        .attr("y", -15)
        .attr("class", "title")
        .text("Comparing Motifs");

}

const prepareDataForVizualization = async (specificDate, dataframe) => {
    let graphData = []
    for (let i = 0; i < dataframe.length; i++) {
        graphData.push({
            index: i,  //  date : specificDate.clone().add(i, 'days')
            value: dataframe[i]
        })
    }

    return graphData;
}

const csvNames = async () => {
    const url = 'http://localhost:3000/csv-files'

    let fileNameContainer = document.getElementById('fileName');
    let fileNameContainerSecond = document.getElementById('dataPath');

    try {
        let fileNames = await axios.get(url);
        for (let fileName of fileNames.data){
            let option = document.createElement("option");
            option.value = fileName.toString().replace('.csv', '') ;
            option.text = fileName.toString().replace('.csv', '');

            fileNameContainer.appendChild(option.cloneNode(true));
            fileNameContainerSecond.appendChild(option.cloneNode(true));
        }

    } catch (error) {
        console.error('File upload failed:', error);
    }
}

const addHtmlSection = async (id) => {
    const targetDiv = document.getElementById("main-container");

    const innerDiv = document.createElement("div");
    innerDiv.classList.add("inner-div");
    innerDiv.id = "id" + id;

    // Create a new div element for the nested div
    const mainChartDiv = document.createElement("div");
    mainChartDiv.classList.add("nested-div");
    mainChartDiv.id = "mainChart" + id;

    const patternChartDiv = document.createElement("div");
    patternChartDiv.classList.add("nested-div");
    patternChartDiv.id = "patternChart" + id;

    const separator = document.createElement("div");
    separator.classList.add("separator");

    // Append the nested div to the inner div
    innerDiv.appendChild(mainChartDiv);
    innerDiv.appendChild(patternChartDiv);
    innerDiv.appendChild(separator);

    targetDiv.appendChild(innerDiv);
}

const createDynamicCheckBox = async (columnToExtractDataContainer, checkBoxName, k) => {
    let label = document.createElement('label');
    label.textContent = (k + 1)
    let checkbox = document.createElement('input');
    checkbox.type = 'checkbox';
    // checkbox.id = k;
    checkbox.value = (k + 1);
    checkbox.name = checkBoxName

    label.appendChild(checkbox)
    columnToExtractDataContainer.append(label)
}

// const showDataSetRow = async (webR, fileName) => {
//     let baseFilePath = 'http://localhost:3000/csv/'
//     let dataa = await webR.evalR(`get.data.bruto("${(baseFilePath + fileName)}")`)
//     let dataframee = await dataa.toJs();
//
//     const headerRow = document.getElementById('headerRow');
//     const tableBody = document.getElementById('tableBody');
//     const columnToExtractDataContainer = document.getElementById('columnToExtractDataId');
//     const columnToFilterContainer = document.getElementById('columnToFilterId');
//
//     let numberOfColumns = dataframee.names.length
//     for (let k = 0; k < numberOfColumns; k++) {
//         createDynamicCheckBox(columnToExtractDataContainer, 'columnToExtractData', k)
//         createDynamicCheckBox(columnToFilterContainer, 'columnToFilter', k)
//     }
//
//
//     for (let columnName of dataframee.names) {
//         let th = document.createElement('th');
//         th.classList.add('table-row');
//         th.textContent = columnName;
//         th.classList.add('table-cell');
//         headerRow.appendChild(th);
//     }
//
//     let numberRows = dataframee.values[0].values.length
//     for (let i = 0; i < numberRows; i++) {
//         let row = document.createElement('row');
//         row.classList.add('table-row');
//         for (let j = 0; j < numberOfColumns; j++) {
//             let td = document.createElement('td');
//             td.textContent = dataframee.values[j].values[i]
//             td.classList.add('table-cell');
//             row.appendChild(td);
//         }
//         tableBody.appendChild(row)
//     }
//
// }

const showDataSetRow = async (webR, fileName) => {
    let baseFilePath = 'http://localhost:3000/csv/'
    let dataa = await webR.evalR(`get.data.bruto("${(baseFilePath + fileName)}")`)
    let dataframee = await dataa.toJs();

    // const headerRow = document.getElementById('headerRow');
    const tableBody = document.getElementById('tableBody');
    const columnToExtractDataContainer = document.getElementById('columnToExtractDataId');
    const columnToFilterContainer = document.getElementById('columnToFilterId');

    let numberOfColumns = dataframee.names.length
    for (let k = 0; k < numberOfColumns; k++) {
        createDynamicCheckBox(columnToExtractDataContainer, 'columnToExtractData', k)
        createDynamicCheckBox(columnToFilterContainer, 'columnToFilter', k)
    }


    let divContainer = document.createElement('div');
    divContainer.classList.add('table-head');
    for (let columnName of dataframee.names) {
        let th = document.createElement('div');
        th.classList.add('table-row');
        th.textContent = columnName;
        th.classList.add('table-cell');
        divContainer.appendChild(th);
    }
    tableBody.appendChild(divContainer);

    let numberRows = dataframee.values[0].values.length
    for (let i = 0; i < numberRows; i++) {
        let row = document.createElement('div');
        row.classList.add('table-row');
        for (let j = 0; j < numberOfColumns; j++) {
            let td = document.createElement('div');
            td.textContent = dataframee.values[j].values[i]
            td.classList.add('table-cell');
            row.appendChild(td);
        }
        tableBody.appendChild(row)
    }

}

// Function to handle the file upload using axios
document.getElementById("uploadFile").onclick = async () => {
    const fileInput = document.getElementById('fileInput');
    const file = fileInput.files[0];

    if (!file) {
        alert('Please select a file to upload.');
        return;
    }

    const formData = new FormData();
    formData.append('file', file);

    // Display the file name in the console (optional)
    console.log('Uploading file:', file.name);
    let url = 'http://localhost:3000/upload/' + file.name;
    let response;
    try {
        response = await axios.post(url, formData, {
            headers: {
                'Content-Type': 'multipart/form-data'
            }
        });

        alert('File uploaded successfully!');
    } catch (error) {
        console.error('File upload failed:', error);
    }
}

const maxArray = async (columnToFilter) => {
    let maxValue = parseInt(columnToFilter[0]);
    for(let mv of columnToFilter){
        if (maxValue <= parseInt(mv)){
            maxValue = parseInt(mv)
        }
    }

    return maxValue;
}

const minArray = async (columnToFilter) => {
    let minValue = parseInt(columnToFilter[0]);
    for(let mv of columnToFilter){
        if (minValue >= parseInt(mv)){
            minValue = parseInt(mv)
        }
    }

    return minValue;
}


document.getElementById("idButton").onclick = async () => {
    let rscript = await exploreMotifScript();

    let dataPath = document.getElementById('dataPath').value;
    let columnToFilter = Array.from(document.querySelectorAll('input[name="columnToFilter"]:checked')).map(checkbox => checkbox.value);
    let filterValue = document.getElementById('filterValue').value;
    let columnToExtractData = Array.from(document.querySelectorAll('input[name="columnToExtractData"]:checked')).map(checkbox => checkbox.value);
    let distanceName = document.getElementById('distanceName').value;
    let sequenceLength = parseInt(document.getElementById('sequenceLength').value);
    let epsilonQuery = parseFloat(document.getElementById('epsilonQuery').value);
    let kk = document.getElementById('k').value;
    let n = parseInt(document.getElementById('n').value);
    let minRow = document.getElementById('min-row').value;
    let maxRow = document.getElementById('max-row').value;

    // You can now send these values to the server or use them in your application logic
    let fileBasePath = 'http://localhost:3000/csv/'
    let getTimeseriesData = `get.timeseries.data("${(fileBasePath + dataPath)}", ${columnToFilter.length == 1 ? columnToFilter[0] : -1}, "${filterValue}", ${columnToExtractData.length > 1 ? await minArray(columnToExtractData) : columnToExtractData[0]}, ${columnToExtractData.length > 1 ? await maxArray(columnToExtractData) : -1}, ${!minRow ? -1 : parseInt(minRow)}, ${!maxRow ? -1 : parseInt(maxRow)})`

    // execute r function
    // let rscript = await exploreMotifScript();
    await webR.evalR(rscript);

    // Loading the time series data
    let data = await webR.evalR(getTimeseriesData)
    let dataframe = await data.toArray();

    // creating time series chart
    // let specificDate = moment('2024-06-08');
    let graphData = await prepareDataForVizualization(null, dataframe);

    let exploreMotif = `explore.motif("${(fileBasePath + dataPath)}", ${columnToFilter.length == 1 ? columnToFilter[0] : -1}, "${filterValue}", ${columnToExtractData.length > 1 ? await minArray(columnToExtractData) : columnToExtractData[0]}, ${sequenceLength}, "${distanceName}", ${epsilonQuery}, ${(parseInt(kk) + 1)}, ${n}, ${columnToExtractData.length > 1 ? await maxArray(columnToExtractData) : -1},  ${!minRow ? -1 : parseInt(minRow)},  ${!maxRow ? -1 : parseInt(maxRow)})`

    // creating similarity matrix and plotting patterns
    let returnValue = await webR.evalR(exploreMotif)
    // let x = await returnValue.toArray();
    let matrix = await returnValue.toJs();
    // one array correspond to one column of matrix
    console.log('n: ', matrix.values[0].values.length)

    let k = parseInt(kk) - 1;
    n = Math.min(matrix.values[0].values.length, n)
    for (let i = k; i < n; i++) {
        let numberElementsNonZero = 0
        let distanceArray = new Array(matrix.values[0].values.length).fill(0);
        for (let j = i; j < matrix.values[i].values.length; j++) {   // let j = k
            if (matrix.values[j].values[i] < epsilonQuery && matrix.values[j].values[i] > 0) {  //
                distanceArray[j] = j
                numberElementsNonZero++
            }
        }
        if (numberElementsNonZero > 0) {
            addHtmlSection(i);
            let mainChartId = "#mainChart" + i;
            let patternChartId = "#patternChart" + i;

            distanceArray = distanceArray.filter(dElement => dElement !== 0)
            grafic(graphData, mainChartId, i, distanceArray, sequenceLength);
            patternVizualzer(graphData, patternChartId, i, distanceArray, sequenceLength)
        }
    }
}

document.getElementById("showDataset").onclick = async () => {
    // loading the script
    let rscript = await exploreMotifScript();
    await webR.evalR(rscript);
    let fileName = document.getElementById('fileName').value;

    showDataSetRow(webR, fileName);
}

csvNames();


