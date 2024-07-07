import * as fs from 'fs';


fs.readFile('r_files/demo.txt', 'utf8', function (err, data) {
    // Display the file content
    console.log(data);
});

console.log('readFile called');