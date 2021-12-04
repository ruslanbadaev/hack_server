// const app = require("express")();
// const http = require("http").Server(app);
// const io = require("socket.io")(http);
// const port = process.env.PORT || 666;

// app.get("/", (req, res) => {
//   res.sendFile(__dirname + "/index.html");
// });

// io.on("connection", (socket) => {
//   socket.on("chat message", (msg) => {
//     console.log(msg);
//     io.emit("chat message", msg);
//   });
// });

// http.listen(port, () => {
//   console.log(`Socket.IO server running at http://localhost:${port}/`);
// });

const http = require("http");
const express = require("express");
const WebSocket = require("ws");
const path = require("path");
const app = express();

const server = http.createServer(app);

const webSocketServer = new WebSocket.Server({ server });
app.use("/", express.static(path.join(__dirname, "dist/")));
app.get("/", (req, res) => {
  res.sendFile(__dirname + "/index.html");
});

app.get("/", (req, res) => {
  res.sendFile(__dirname + "/dist/bundle.main.js");
});
app.get("/", (req, res) => {
  res.sendFile(__dirname + "/bundle.main.js");
});

webSocketServer.on("connection", (ws) => {
  ws.on("message", (m) => {
    console.log(JSON.parse(m.toString()));
    ws.send(m.buffer);
    webSocketServer.clients.forEach((client) => {
      client.send(m);
      // ws.send(m);
    });
  });

  ws.on("error", (e) => ws.send(e));

  ws.send("Hi there, I am a WebSocket server");
});

server.listen(3000, () => console.log("Server started"));
