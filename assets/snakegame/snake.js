(() => {
  const canvas = document.getElementById("game");
  const ctx = canvas.getContext("2d");
  const pauseBtn = document.getElementById("pause");

  // Hücre boyutu (px)
  const size = 24;
  let cellsX, cellsY;
  let snake, dir, nextDir, food, score, best, speed, tick, paused, dead;

  function resizeCanvas() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
    cellsX = Math.max(8, Math.floor(canvas.width / size));
    cellsY = Math.max(8, Math.floor(canvas.height / size));
  }
  window.addEventListener("resize", resizeCanvas, { passive: true });
  resizeCanvas();

  const rndCellX = () => Math.floor(Math.random() * cellsX);
  const rndCellY = () => Math.floor(Math.random() * cellsY);

  function reset() {
    snake = [{ x: Math.floor(cellsX / 2), y: Math.floor(cellsY / 2) }];
    dir = { x: 1, y: 0 };
    nextDir = { x: 1, y: 0 };
    spawnFood();
    score = 0;
    speed = 7;
    tick = 0;
    paused = false;
    dead = false;
    drawScore();
    drawFrame();
  }

  function spawnFood() {
    do { food = { x: rndCellX(), y: rndCellY() }; }
    while (snake.some(s => s.x === food.x && s.y === food.y));
  }

  function drawGrid() {
    ctx.fillStyle = "#0f1720";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
  }

  function drawSnake() {
    snake.forEach((s, i) => {
      ctx.fillStyle = i === 0 ? "#94e2ff" : "#4cc9f0";
      ctx.fillRect(s.x * size, s.y * size, size, size);
    });
  }

  function drawFood() {
    ctx.fillStyle = "#f97316";
    ctx.fillRect(food.x * size, food.y * size, size, size);
  }

  function drawFrame() {
    drawGrid(); drawFood(); drawSnake();
  }

  function setDirectionByName(name) {
    const map = {
      up:    { x: 0,  y: -1 },
      down:  { x: 0,  y: 1 },
      left:  { x: -1, y: 0 },
      right: { x: 1,  y: 0 }
    };
    const nd = map[name];
    if (!nd) return;
    // Anında 180° dönmeyi engelle
    if (nd.x !== -dir.x || nd.y !== -dir.y) nextDir = nd;
  }

  function step() {
    if (paused || dead) return;
    requestAnimationFrame(step);

    // hız kontrol
    if (++tick < Math.max(1, Math.floor(60 / speed))) return;
    tick = 0;

    if ((nextDir.x !== -dir.x) || (nextDir.y !== -dir.y)) dir = nextDir;

    const head = { x: snake[0].x + dir.x, y: snake[0].y + dir.y };

    // wrap-around (duvara çarpınca diğer taraftan çık)
    if (head.x < 0) head.x = cellsX - 1;
    if (head.x >= cellsX) head.x = 0;
    if (head.y < 0) head.y = cellsY - 1;
    if (head.y >= cellsY) head.y = 0;

    // kendine çarpma → oyun biter
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

    drawFrame();
  }

  function drawScore() {
    document.getElementById("score").textContent = `Score: ${score}`;
    best = Math.max(score, Number(localStorage.getItem("snake_best") || 0));
    localStorage.setItem("snake_best", String(best));
    document.getElementById("best").textContent = `Best: ${best}`;
  }

  function gameOver() {
    dead = true;
    drawFrame();
    ctx.fillStyle = "rgba(0,0,0,0.55)";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = "#e8eef6";
    ctx.font = "bold 32px system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Game Over", canvas.width / 2, canvas.height / 2);
  }

  // --- Klavye ---
  const keyMap = {
    ArrowUp: { x: 0, y: -1 }, KeyW: { x: 0, y: -1 },
    ArrowDown: { x: 0, y: 1 }, KeyS: { x: 0, y: 1 },
    ArrowLeft: { x: -1, y: 0 }, KeyA: { x: -1, y: 0 },
    ArrowRight: { x: 1, y: 0 }, KeyD: { x: 1, y: 0 },
  };
  addEventListener("keydown", (e) => {
    if (e.code === "KeyP") { togglePause(); return; }
    if (dead) return;
    const nd = keyMap[e.code];
    if (nd && (nd.x !== -dir.x || nd.y !== -dir.y)) nextDir = nd;
  });

  // --- D-Pad: tık ve dokunma ---
  document.querySelectorAll(".btn-dir").forEach(btn => {
    const handler = (ev) => {
      ev.preventDefault();
      if (dead) return;
      setDirectionByName(btn.dataset.dir);
    };
    btn.addEventListener("click", handler);
    btn.addEventListener("touchstart", handler, { passive: false });
  });

  // --- Swipe (kaydırma) desteği ---
  let touchStartX = 0, touchStartY = 0, touchStartT = 0;
  const SWIPE_MIN = 24; // px

  canvas.addEventListener("touchstart", (e) => {
    if (!e.touches || e.touches.length === 0) return;
    const t = e.touches[0];
    touchStartX = t.clientX; touchStartY = t.clientY; touchStartT = performance.now();
  }, { passive: true });

  canvas.addEventListener("touchend", (e) => {
    const t = (e.changedTouches && e.changedTouches[0]) ? e.changedTouches[0] : null;
    if (!t) return;
    const dx = t.clientX - touchStartX;
    const dy = t.clientY - touchStartY;
    if (Math.abs(dx) < SWIPE_MIN && Math.abs(dy) < SWIPE_MIN) return;

    if (Math.abs(dx) > Math.abs(dy)) {
      setDirectionByName(dx > 0 ? "right" : "left");
    } else {
      setDirectionByName(dy > 0 ? "down" : "up");
    }
  }, { passive: true });

  // İsteğe bağlı: mouse sürükleme ile yön (PC)
  let mouseDown = false, mx0 = 0, my0 = 0;
  canvas.addEventListener("mousedown", (e) => { mouseDown = true; mx0 = e.clientX; my0 = e.clientY; });
  canvas.addEventListener("mouseup", (e) => {
    if (!mouseDown) return;
    mouseDown = false;
    const dx = e.clientX - mx0, dy = e.clientY - my0;
    if (Math.abs(dx) < SWIPE_MIN && Math.abs(dy) < SWIPE_MIN) return;
    if (Math.abs(dx) > Math.abs(dy)) {
      setDirectionByName(dx > 0 ? "right" : "left");
    } else {
      setDirectionByName(dy > 0 ? "down" : "up");
    }
  });

  // Pause/Play
  function togglePause() {
    paused = !paused;
    if (!paused && !dead) requestAnimationFrame(step);
  }
  pauseBtn.addEventListener("click", (e) => { e.preventDefault(); togglePause(); });
  pauseBtn.addEventListener("touchstart", (e) => { e.preventDefault(); togglePause(); }, { passive: false });

  // Restart
  document.getElementById("restart").addEventListener("click", (e) => { e.preventDefault(); reset(); requestAnimationFrame(step); });
  document.getElementById("restart").addEventListener("touchstart", (e) => { e.preventDefault(); reset(); requestAnimationFrame(step); }, { passive: false });

  // Başlat
  reset();
  requestAnimationFrame(step);
})();
