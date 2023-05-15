function lin(n: number) {
    return Array(n).fill(0).map((x, i) => i)
}

function logspace(x: number, y: number, len: number) {
    var arr, end, tmp, d;
    var a = Math.log10(x)
    var b = Math.log10(y)

    if (len <= 0) {
        return [];
    }

    end = len - 1;
    d = (b - a) / end;

    arr = new Array(len);
    tmp = a;
    arr[0] = Math.pow(10, tmp);
    for (var i = 1; i < end; i++) {
        tmp += d;
        arr[i] = Math.pow(10, tmp);
    }
    arr[end] = Math.pow(10, b);
    return arr;
}

function lnmode(A: number[], x: number[]) {
    let prefactor = A[0] / (Math.sqrt(2 * Math.PI) * Math.log(A[2]))
    let div = 2 * Math.pow(Math.log(A[2]), 2)
    let psd = x.map(x => prefactor * Math.exp(-Math.pow(Math.log(x / A[1]), 2) / div))
    return psd
}