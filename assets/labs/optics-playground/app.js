/* Optics Playground — all client-side, no deps. */
(() => {
  /* ----------- Utils ----------- */
  const $ = sel => document.querySelector(sel);
  const $$ = sel => document.querySelectorAll(sel);
  const clamp = (v, a, b) => Math.max(a, Math.min(b, v));

  function setCanvasSize(c, w, h) {
    const dpr = Math.max(1, Math.min(2, window.devicePixelRatio || 1));
    c.width = Math.round(w * dpr);
    c.height = Math.round(h * dpr);
    c.style.width = `${w}px`;
    c.style.height = `${h}px`;
    const ctx = c.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return ctx;
  }

  function drawAxes(ctx, w, h, label) {
    ctx.clearRect(0,0,w,h);
    ctx.fillStyle = '#0a121a';
    ctx.fillRect(0,0,w,h);
    ctx.strokeStyle = '#223242';
    ctx.lineWidth = 1;
    // frame
    ctx.strokeRect(0.5,0.5,w-1,h-1);
    // ticks
    ctx.fillStyle = '#9fb0c3';
    ctx.font = '12px system-ui';
    ctx.textAlign = 'left';
    ctx.fillText(label || '', 6, 14);
  }

  /* ----------- 1) GSD Explorer ----------- */
  const alt = $('#alt_km'), fmm = $('#f_mm'), pix = $('#pix_um'), nx = $('#nx'), ny = $('#ny');
  const gsdOut = $('#gsd_out'), swx = $('#swath_x'), swy = $('#swath_y'), fvx = $('#fov_x'), fvy = $('#fov_y');
  const gsdC = $('#gsdPlot'); const gsdCtx = gsdC.getContext('2d');

  function computeGSD() {
    const H = Number(alt.value) * 1e3;         // m
    const f = Number(fmm.value) * 1e-3;        // m
    const p = Number(pix.value) * 1e-6;        // m
    const NX = Number(nx.value), NY = Number(ny.value);

    const gsd = H * (p / f);                   // m/px
    const sensX = NX * p, sensY = NY * p;      // m
    const fovX = 2 * Math.atan((sensX/2) / f); // rad
    const fovY = 2 * Math.atan((sensY/2) / f); // rad
    const swathX = 2 * H * Math.tan(fovX/2);   // m
    const swathY = 2 * H * Math.tan(fovY/2);   // m

    gsdOut.textContent = gsd.toFixed(3);
    swx.textContent = (swathX/1000).toFixed(2);
    swy.textContent = (swathY/1000).toFixed(2);
    fvx.textContent = (fovX*180/Math.PI).toFixed(2);
    fvy.textContent = (fovY*180/Math.PI).toFixed(2);

    // footprint preview (not to scale — normalized)
    const w = gsdC.clientWidth || 480, h = gsdC.clientHeight || 240;
    setCanvasSize(gsdC, w, h);
    drawAxes(gsdCtx, w, h, 'Footprint');
    const pad = 20, maxRect = Math.min(w, h) - pad*2;
    const ratio = (swathX / swathY) || 1;
    let rw = maxRect, rh = maxRect;
    if (ratio > 1) rh = maxRect / ratio; else rw = maxRect * ratio;
    gsdCtx.fillStyle = '#4cc9f0';
    gsdCtx.globalAlpha = 0.15;
    gsdCtx.fillRect((w-rw)/2, (h-rh)/2, rw, rh);
    gsdCtx.globalAlpha = 1;
    gsdCtx.strokeStyle = '#4cc9f0';
    gsdCtx.lineWidth = 2;
    gsdCtx.strokeRect((w-rw)/2, (h-rh)/2, rw, rh);
    gsdCtx.fillStyle = '#9fb0c3';
    gsdCtx.font = '12px system-ui';
    gsdCtx.textAlign = 'center';
    gsdCtx.fillText(`~ ${ (swathX/1000).toFixed(1) } km × ${ (swathY/1000).toFixed(1) } km`, w/2, (h+rh)/2 + 16);
  }

  $('#btn-gsd').addEventListener('click', computeGSD);
  $('#auto-gsd').addEventListener('change', computeGSD);
  [alt,fmm,pix,nx,ny].forEach(el => el.addEventListener('input', () => { if ($('#auto-gsd').checked) computeGSD(); }));
  computeGSD();

  /* ----------- FFT (1D & 2D) ----------- */
  // In-place radix-2 Cooley–Tukey. re[], im[] are Float64Array. n is power-of-two.
  function fft1d(re, im, inverse=false) {
    const n = re.length;
    // bit-reversal
    for (let i=0, j=0; i<n; i++) {
      if (i<j) { let tr=re[i]; re[i]=re[j]; re[j]=tr; tr=im[i]; im[i]=im[j]; im[j]=tr; }
      let m = n>>1; while (m>=1 && j>=m) { j-=m; m>>=1; } j+=m;
    }
    for (let len=2; len<=n; len<<=1) {
      const ang = 2*Math.PI/len * (inverse ? -1 : 1);
      const wlen_r = Math.cos(ang), wlen_i = Math.sin(ang);
      for (let i=0; i<n; i+=len) {
        let wr=1, wi=0;
        for (let j=0; j<(len>>1); j++) {
          const u_r = re[i+j], u_i = im[i+j];
          const v_r = re[i+j+(len>>1)]*wr - im[i+j+(len>>1)]*wi;
          const v_i = re[i+j+(len>>1)]*wi + im[i+j+(len>>1)]*wr;
          re[i+j] = u_r + v_r;  im[i+j] = u_i + v_i;
          re[i+j+(len>>1)] = u_r - v_r;  im[i+j+(len>>1)] = u_i - v_i;
          const nxt_wr = wr*wlen_r - wi*wlen_i;
          wi = wr*wlen_i + wi*wlen_r; wr = nxt_wr;
        }
      }
    }
    if (inverse) for (let i=0;i<n;i++){ re[i]/=n; im[i]/=n; }
  }

  function fft2d(re, im, N, inverse=false) {
    // rows
    let rowRe = new Float64Array(N), rowIm = new Float64Array(N);
    for (let y=0;y<N;y++){
      for (let x=0;x<N;x++){ rowRe[x]=re[y*N+x]; rowIm[x]=im[y*N+x]; }
      fft1d(rowRe, rowIm, inverse);
      for (let x=0;x<N;x++){ re[y*N+x]=rowRe[x]; im[y*N+x]=rowIm[x]; }
    }
    // cols
    let colRe = new Float64Array(N), colIm = new Float64Array(N);
    for (let x=0;x<N;x++){
      for (let y=0;y<N;y++){ colRe[y]=re[y*N+x]; colIm[y]=im[y*N+x]; }
      fft1d(colRe, colIm, inverse);
      for (let y=0;y<N;y++){ re[y*N+x]=colRe[y]; im[y*N+x]=colIm[y]; }
    }
  }

  function fftshift2d(arr, N) {
    // quadrant swap
    const half = N>>1;
    for (let y=0;y<half;y++){
      for (let x=0;x<half;x++){
        const a = y*N + x;
        const b = (y+half)*N + (x+half);
        const c = y*N + (x+half);
        const d = (y+half)*N + x;
        let t = arr[a]; arr[a] = arr[b]; arr[b] = t;
        t = arr[c]; arr[c] = arr[d]; arr[d] = t;
      }
    }
  }

  /* ----------- 2) Diffraction & PSF/MTF ----------- */
  const Nsel = $('#N'), eps = $('#eps'), arms = $('#arms'), wsp = $('#wsp'), ang = $('#ang');
  const eps_v = $('#eps_v'), arms_v = $('#arms_v'), wsp_v = $('#wsp_v'), ang_v = $('#ang_v');
  const psfCanvas = $('#psf'), mtfCanvas = $('#mtf');

  [eps,arms,wsp,ang].forEach(el => el.addEventListener('input', () => {
    eps_v.textContent = Number(eps.value).toFixed(2);
    arms_v.textContent = arms.value;
    wsp_v.textContent = `${Number(wsp.value).toFixed(1)}%`;
    ang_v.textContent = `${ang.value}°`;
    if ($('#auto-psf').checked) renderPSF();
  }));
  Nsel.addEventListener('change', () => { if ($('#auto-psf').checked) renderPSF(); });
  $('#btn-psf').addEventListener('click', renderPSF);
  $('#auto-psf').addEventListener('change', renderPSF);

  function buildPupil(N, epsilon, arms, spiderPct, rotDeg) {
    const re = new Float64Array(N*N);
    const im = new Float64Array(N*N);
    const R = (N/2 - 2); // aperture radius in px
    const epsR = epsilon * R;
    const rot = rotDeg * Math.PI/180;
    const spiderW = (spiderPct/100) * R; // thickness (px)
    const cx = (N-1)/2, cy = (N-1)/2;

    for (let y=0;y<N;y++){
      const yy = y - cy;
      for (let x=0;x<N;x++){
        const xx = x - cx;
        const r = Math.hypot(xx, yy);
        let amp = (r<=R && r>=epsR) ? 1 : 0;

        // spiders: arms equally spaced
        if (amp && arms>0 && spiderW>0){
          for (let k=0;k<arms;k++){
            const phi = rot + k * Math.PI/arms;
            const dist = Math.abs( xx*Math.sin(phi) - yy*Math.cos(phi) ); // distance to line
            if (dist <= spiderW*0.5) { amp = 0; break; }
          }
        }

        // pre-center (multiply by (-1)^(x+y)) to put DC at center after FFT
        const sign = ((x+y)&1) ? -1 : 1;
        re[y*N+x] = amp * sign;
        im[y*N+x] = 0;
      }
    }
    return {re, im};
  }

  function renderPSF() {
    const N = Number(Nsel.value);
    const epsilon = Number(eps.value);
    const armsN = Number(arms.value);
    const spiderPct = Number(wsp.value);
    const rotDeg = Number(ang.value);

    // pupil
    let {re, im} = buildPupil(N, epsilon, armsN, spiderPct, rotDeg);

    // field = FFT(pupil)
    fft2d(re, im, N, false);

    // intensity
    const I = new Float64Array(N*N);
    let maxI = 0;
    for (let i=0;i<I.length;i++){ const v = re[i]*re[i] + im[i]*im[i]; I[i]=v; if (v>maxI) maxI=v; }

    // draw PSF (log)
    drawImageFromArray(psfCanvas, I, N, true);

    // OTF = FFT(PSF)  -> magnitude -> radial MTF
    let R2 = I.slice(); // copy
    let Z = new Float64Array(R2.length);
    // pre-center PSF too
    for (let y=0;y<N;y++) for (let x=0;x<N;x++){ const idx=y*N+x; const s=((x+y)&1)?-1:1; R2[idx]*=s; }
    fft2d(R2, Z, N, false);
    for (let i=0;i<R2.length;i++){ R2[i] = Math.hypot(R2[i], Z[i]); }
    // normalize by DC
    const DC = R2[((N>>1))*N + (N>>1)] || 1;
    for (let i=0;i<R2.length;i++){ R2[i] /= DC; }

    // radial average to 1D curve
    const mtfCurve = radialAverage(R2, N);
    drawMTF(mtfCanvas, mtfCurve);
  }

  function drawImageFromArray(canvas, arr, N, logScale=false) {
    const w = canvas.clientWidth || 256, h = canvas.clientHeight || 256;
    const ctx = setCanvasSize(canvas, w, h);
    // normalize
    let maxV = 0;
    for (let i=0;i<arr.length;i++) if (arr[i]>maxV) maxV=arr[i];
    const img = ctx.createImageData(N, N);
    for (let i=0;i<arr.length;i++){
      let v = arr[i]/(maxV||1);
      if (logScale) v = Math.log10(1 + 1000*v)/Math.log10(1001); // gentle log
      const g = Math.round(255*clamp(v,0,1));
      const j = i*4; img.data[j]=img.data[j+1]=img.data[j+2]=g; img.data[j+3]=255;
    }
    // nearest-neighbor upscale
    const off = new OffscreenCanvas(N, N);
    const octx = off.getContext('2d');
    octx.putImageData(img,0,0);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(off, 0,0, w, h);
  }

  function radialAverage(arr, N) {
    const cx = (N-1)/2, cy = (N-1)/2;
    const maxR = Math.floor(Math.hypot(cx, cy));
    const bins = new Float64Array(maxR+1);
    const cnts = new Uint32Array(maxR+1);
    for (let y=0;y<N;y++){
      for (let x=0;x<N;x++){
        const r = Math.round(Math.hypot(x-cx, y-cy));
        if (r<=maxR){ bins[r]+=arr[y*N+x]; cnts[r]++; }
      }
    }
    const curve = new Float64Array(maxR+1);
    for (let r=0;r<=maxR;r++){ curve[r] = cnts[r] ? (bins[r]/cnts[r]) : 0; }
    // normalize to 1 at r=0
    const dc = curve[0]||1;
    for (let r=0;r<curve.length;r++) curve[r]/=dc;
    return curve;
  }

  function drawMTF(canvas, curve) {
    const w = canvas.clientWidth || 256, h = canvas.clientHeight || 256;
    const ctx = setCanvasSize(canvas, w, h);
    drawAxes(ctx, w, h, 'MTF (normalized spatial freq)');
    // plot
    ctx.strokeStyle = '#4cc9f0';
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i=0;i<curve.length;i++){
      const x = (i/(curve.length-1)) * (w-16) + 8;
      const y = (1-curve[i]) * (h-16) + 8;
      if (i===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
    }
    ctx.stroke();
  }

  renderPSF();

  /* ----------- 3) Misalignment via Zernike ----------- */
  const NZ = $('#NZ'), epsZ = $('#epsZ');
  const psfZ = $('#psfZ'), phaseZ = $('#phaseZ');
  const zCtrls = ['z_tiltx','z_tilty','z_def','z_ast0','z_ast45','z_comax','z_comay'];
  zCtrls.forEach(id => {
    const el = $('#'+id), lab = $('#'+id+'_v');
    el.addEventListener('input', () => { lab.textContent = Number(el.value).toFixed(2); if ($('#auto-z').checked) renderZ(); });
    lab.textContent = Number(el.value).toFixed(2);
  });
  epsZ.addEventListener('input', () => { $('#epsZ_v').textContent = Number(epsZ.value).toFixed(2); if ($('#auto-z').checked) renderZ(); });
  NZ.addEventListener('change', () => { if ($('#auto-z').checked) renderZ(); });
  $('#btn-z').addEventListener('click', renderZ);
  $('#auto-z').addEventListener('change', renderZ);

  function buildPupilZ(N, epsilon, coeffs) {
    const re = new Float64Array(N*N);
    const im = new Float64Array(N*N);
    const phase = new Float64Array(N*N);
    const R = (N/2 - 2);
    const epsR = epsilon * R;
    const cx = (N-1)/2, cy = (N-1)/2;

    for (let y=0;y<N;y++){
      const yy = y - cy;
      for (let x=0;x<N;x++){
        const xx = x - cx;
        const rpx = Math.hypot(xx, yy);
        const idx = y*N+x;
        if (rpx>R || rpx<epsR) { re[idx]=im[idx]=0; phase[idx]=0; continue; }

        // normalized polar coords in unit radius
        const r = (rpx - epsR) / (R - epsR); // map annulus to [0,1]
        const th = Math.atan2(yy, xx);

        // Zernikes (approx normalization; coeff in waves)
        const Z_tiltx = 2 * r * Math.cos(th);
        const Z_tilty = 2 * r * Math.sin(th);
        const Z_def   = Math.sqrt(3) * (2*r*r - 1);
        const Z_ast0  = Math.sqrt(6) * r*r * Math.cos(2*th);
        const Z_ast45 = Math.sqrt(6) * r*r * Math.sin(2*th);
        const Z_comax = Math.sqrt(8) * (3*r*r*r - 2*r) * Math.cos(th);
        const Z_comay = Math.sqrt(8) * (3*r*r*r - 2*r) * Math.sin(th);

        const waves =
          coeffs.tiltx*Z_tiltx + coeffs.tilty*Z_tilty + coeffs.def*Z_def +
          coeffs.ast0*Z_ast0 + coeffs.ast45*Z_ast45 + coeffs.comax*Z_comax + coeffs.comay*Z_comay;

        const ph = 2*Math.PI * waves;
        phase[idx] = ph;

        const sign = ((x+y)&1) ? -1 : 1;
        re[idx] = Math.cos(ph) * sign;
        im[idx] = Math.sin(ph) * sign;
      }
    }
    return {re, im, phase};
  }

  function renderZ() {
    const N = Number(NZ.value);
    const epsilon = Number(epsZ.value);
    const coeffs = {
      tiltx: Number($('#z_tiltx').value),
      tilty: Number($('#z_tilty').value),
      def:   Number($('#z_def').value),
      ast0:  Number($('#z_ast0').value),
      ast45: Number($('#z_ast45').value),
      comax: Number($('#z_comax').value),
      comay: Number($('#z_comay').value),
    };

    let {re, im, phase} = buildPupilZ(N, epsilon, coeffs);

    // field = FFT(pupil with phase)
    fft2d(re, im, N, false);

    // intensity
    const I = new Float64Array(N*N);
    for (let i=0;i<I.length;i++) I[i] = re[i]*re[i] + im[i]*im[i];

    drawImageFromArray(psfZ, I, N, true);
    drawPhase(phaseZ, phase, N);
  }

  function drawPhase(canvas, phase, N) {
    const w = canvas.clientWidth || 256, h = canvas.clientHeight || 256;
    const ctx = setCanvasSize(canvas, w, h);
    // map phase [-π, π] to gradient
    const img = ctx.createImageData(N, N);
    for (let i=0;i<phase.length;i++){
      let ph = phase[i];
      while (ph > Math.PI) ph -= 2*Math.PI;
      while (ph < -Math.PI) ph += 2*Math.PI;
      const v = (ph + Math.PI) / (2*Math.PI); // 0..1
      // simple blue→orange map
      const r = Math.round(255 * v);
      const g = Math.round(255 * (0.5 + 0.2*Math.sin(6.28*(v-0.25))));
      const b = Math.round(255 * (1-v));
      const j = i*4; img.data[j]=r; img.data[j+1]=g; img.data[j+2]=b; img.data[j+3]=255;
    }
    const off = new OffscreenCanvas(N, N);
    const octx = off.getContext('2d');
    octx.putImageData(img,0,0);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(off, 0,0, w, h);
  }

  renderZ();

  /* Kick responsive redraw on resize */
  window.addEventListener('resize', () => {
    computeGSD();
    renderPSF();
    renderZ();
  }, { passive: true });
})();
