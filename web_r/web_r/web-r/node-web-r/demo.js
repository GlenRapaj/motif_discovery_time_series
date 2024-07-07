import { WebR } from 'webr';
import * as fs from 'fs';
import { Canvas } from 'canvas';


const webR = new WebR();
await webR.init();

let result = await webR.evalR('rnorm(10,5,1)');
let output = await result.toArray();

console.log('Result of running `rnorm` from webR: ', output);
console.log('Result of running `rnorm` from webR result : ', result);


// executing r function in js

let sin = await webR.evalR('sin');
let resultt = await sin([1, 2, 3]);

console.log('Result of running `await sin([1,2,3])`:')
console.log(resultt.values);

