First install Rust if not already done:
```bash
curl https://sh.rustup.rs -sSf | sh
```
(More details [here](https://www.rust-lang.org/tools/install) if you are not using Linux.)

To build the WASM and Javascript files, you need to install [wasm-bindgen](https://github.com/rustwasm/wasm-bindgen) in addition to the `nightly` version of Rust:

```bash
rustup target add wasm32-unknown-unknown --toolchain nightly
cargo +nightly install wasm-bindgen-cli
```

To test the demo page, you have to load the WASM/Javascript files from an HTTP server.
To start a local HTTP server, you can use (requires python):

```bash
# For python 2
python -m SimpleHTTPServer
# For python 3
python -m http.server
```

Then load in your favourite browser the URL [http://localhost:8000/test.html](http://localhost:8000/test.html) and/or [http://localhost:8000/demo.html](http://localhost:8000/demo.html)

