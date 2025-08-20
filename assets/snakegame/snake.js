(() => {
  const canvas = document.getElementById("game");
  const ctx = canvas.getContext("2d");

  const size = 24;                 // grid kare boyutu
  const cells = canvas.width / size; // 480/24 = 20x20
  let snake, dir, nextDir, food, score, best, speed, tick, paused, dead;

  const rndCell = () => Math.floor(Math.random() * cells);

  function reset() {
    snake = [{ x: 10, y: 10 }];
    dir = { x: 1, y: 0 };
    nextDir = { x: 1, y: 0 };
    spawnFood();
    score = 0;
    speed = 7;          // başlangıç hızı
    tick = 0;
    paused = false;
    dead = false;
    drawScore();
  }

  function spawnFood() {
    do { food = { x: rndCell(), y: rndCell() }; }
    while (snake.some(s => s.x === food.x && s.y === food.y));
  }

  function drawGrid() {
    ctx.fillStyle = "#0f1720";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.strokeStyle = "rgba(255,255,255,0.05)";
    for (let i = 1; i < cells; i++) {
      ctx.beginPath(); ctx.moveTo(i * size, 0); ctx.lineTo(i * size, canvas.height); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(0, i * size); ctx.lineTo(canvas.width, i * size); ctx.stroke();
    }
  }

  function drawSnake() {
    snake.forEach((s, i) => {
      const isHead = i === 0;
      ctx.fillStyle = isHead ? "#94e2ff" : "#4cc9f0";
      ctx.fillRect(s.x * size + 2, s.y * size + 2, size - 4, size - 4);
    });
  }

  function drawFood() {
    ctx.fillStyle = "#f97316";
    const pad = 4;
    ctx.fillRect(food.x * size + pad, food.y * size + pad, size - pad * 2, size - pad * 2);
  }

  function step() {
    if (paused || dead) return;
    requestAnimationFrame(step);

    // hız kontrol: belirli aralıklarla güncelle
    if (++tick < Math.max(1, Math.floor(60 / speed))) return;
    tick = 0;

    // yön güncelle
    if ((nextDir.x !== -dir.x) || (nextDir.y !== -dir.y)) dir = nextDir;

    const head = { x: snake[0].x + dir.x, y: snake[0].y + dir.y };

    // duvar çarpışması
    if (head.x < 0 || head.y < 0 || head.x >= cells || head.y >= cells) return gameOver();

    // kendine çarpma
    if (snake.some(s => s.x === head.x && s.y === head.y)) return gameOver();

    snake.unshift(head);

    if (head.x === food.x && head.y === food.y) {
      score++;
      speed = Math.min(20, speed + 0.3);
      spawnFood();
      drawScore();
    } else {
      snake.pop();
    }

    drawGrid(); drawFood(); drawSnake();
  }

  function drawScore() {
    document.getElementById("score").textContent = `Score: ${score}`;
    best = Math.max(score, Number(localStorage.getItem("snake_best") || 0));
    localStorage.setItem("snake_best", String(best));
    document.getElementById("best").textContent = `Best: ${best}`;
  }

  function gameOver() {
    dead = true;
    drawGrid(); drawFood(); drawSnake();
    ctx.fillStyle = "rgba(0,0,0,0.55)";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = "#e8eef6";
    ctx.font = "bold 24px system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Game Over", canvas.width / 2, canvas.height / 2 - 10);
    ctx.font = "16px system-ui, sans-serif";
    ctx.fillText("Restart için butona bas", canvas.width / 2, canvas.height / 2 + 16);
  }

  // Klavye
  const keyMap = {
    ArrowUp: { x: 0, y: -1 }, KeyW: { x: 0, y: -1 },
    ArrowDown: { x: 0, y: 1 }, KeyS: { x: 0, y: 1 },
    ArrowLeft: { x: -1, y: 0 }, KeyA: { x: -1, y: 0 },
    ArrowRight: { x: 1, y: 0 }, KeyD: { x: 1, y: 0 },
  };
  addEventListener("keydown", (e) => {
    if (e.code === "KeyP") { paused = !paused; if (!paused && !dead) requestAnimationFrame(step); return; }
    if (dead) return;
    const nd = keyMap[e.code];
    if (nd) nextDir = nd;
  });

  // Restart
  document.getElementById("restart").addEventListener("click", () => { reset(); requestAnimationFrame(step); });

  // Başlat
  reset();
  drawGrid(); drawFood(); drawSnake();
  requestAnimationFrame(step);
})();

