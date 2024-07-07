import * as fs from 'fs';
import express from 'express';
import cors from 'cors';
import * as path from 'path';
import * as util from 'util';

const app = express()
const port = 3000
app.use(cors())

const dir = process.cwd()
app.use(express.static(path.join(dir, 'public')));


app.get('/', (req, res) => {
    res.sendFile(path.join(dir, 'public', 'index.html'));
});

app.get('/home', (req, res) => {
    res.sendFile(path.join(dir, 'public', 'home.html'));
});

app.get('/about', (req, res) => {
    res.sendFile(path.join(dir, 'public', 'about.html'));
});

app.post('/upload/:fileName', async (req, res) => {
    let data = [];
    req.on("data", (chunk) => {
        data.push(chunk);
    });

    let fileName = req.params.fileName;

    let dir = process.cwd()
    const filePath = `${dir}/r_files/${fileName}`;

    req.on("end", () => {
        let fileData = Buffer.concat(data);
        fs.writeFile(filePath, fileData, "base64", (err) => {
            if (err) {
                res.statusCode = 500;
            }
        });
    });
});


// app.get('/csv',async  (req, res) => {
//     // Advertising_Data.csv
//     let dir = process.cwd()
//     const file = `${dir}/r_files/Advertising_Data.csv`;
//     // res.send('Hello World!')
//     res.download(file);
// })

app.get('/csv/:fileName', async (req, res) => {
    let dir = process.cwd()
    const file = `${dir}/r_files/${req.params.fileName}.csv`;
    res.download(file);
})

app.get('/csv-files', async (req, res) => {
    const readdirPromise = util.promisify(fs.readdir);

    let dir = process.cwd()
    const directoryPath = `${dir}/r_files/`;
    const fileNames = await readdirPromise(directoryPath);

    res.json(fileNames);

})

app.listen(port, () => {
    console.log(`Example app listening on port ${port}`)
})

