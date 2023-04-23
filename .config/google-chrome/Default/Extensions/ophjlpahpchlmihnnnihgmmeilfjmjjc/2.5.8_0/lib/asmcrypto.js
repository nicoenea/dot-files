!function(t, e) {
    "use strict";
    function r() {
        var t = Error.apply(this, arguments);
        this.message = t.message, this.stack = t.stack;
    }
    function i() {
        var t = Error.apply(this, arguments);
        this.message = t.message, this.stack = t.stack;
    }
    function n() {
        var t = Error.apply(this, arguments);
        this.message = t.message, this.stack = t.stack;
    }
    function s(t) {
        for (var e = t.length, r = new Uint8Array(e), i = 0; e > i; i++) {
            var n = t.charCodeAt(i);
            if (n >>> 8) throw new Error("Wide characters are not allowed");
            r[i] = n;
        }
        return r;
    }
    function a(t) {
        for (var e = "", r = 0; r < t.length; r++) e += String.fromCharCode(t[r]);
        return e;
    }
    function h(t) {
        for (var e = "", r = 0; r < t.length; r++) {
            var i = (255 & t[r]).toString(16);
            i.length < 2 && (e += "0"), e += i;
        }
        return e;
    }
    function o(t) {
        return btoa(a(t));
    }
    function f(t) {
        return t -= 1, t |= t >>> 1, t |= t >>> 2, t |= t >>> 4, t |= t >>> 8, (t |= t >>> 16) + 1;
    }
    function u(t) {
        return "number" == typeof t;
    }
    function l(t) {
        return "string" == typeof t;
    }
    function c(t) {
        return t instanceof ArrayBuffer;
    }
    function w(t) {
        return t instanceof Uint8Array;
    }
    function y(t) {
        return t instanceof Int8Array || t instanceof Uint8Array || t instanceof Int16Array || t instanceof Uint16Array || t instanceof Int32Array || t instanceof Uint32Array || t instanceof Float32Array || t instanceof Float64Array;
    }
    function p(t, e, r) {
        "use asm";
        var i = 0, n = 0, s = 0, a = 0, h = 0, o = 0, f = 0, u = 0, l = 0, c = 0, w = 0, y = 0, p = 0, v = 0, d = 0, g = 0;
        var b = 0;
        var m = 0, S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0;
        var L = 0, H = 0, I = 0, q = 0, M = 0, C = 0, Z = 0, z = 0, O = 0, T = 0, R = 0, B = 0, P = 0, K = 0, N = 0, j = 0, D = 0, G = 0, F = 0, W = 0, V = 0, J = 0, Q = 0, X = 0, Y = 0, $ = 0, tt = 0, et = 0, rt = 0, it = 0, nt = 0, st = 0, at = 0, ht = 0, ot = 0, ft = 0, ut = 0, lt = 0, ct = 0, wt = 0, yt = 0, pt = 0, vt = 0, dt = 0, gt = 0, bt = 0, mt = 0, St = 0, At = 0, Et = 0, _t = 0, kt = 0, Ut = 0, xt = 0, Lt = 0, Ht = 0, It = 0, qt = 0, Mt = 0, Ct = 0, Zt = 0, zt = 0, Ot = 0, Tt = 0, Rt = 0, Bt = 0, Pt = 0, Kt = 0, Nt = 0, jt = 0, Dt = 0, Gt = 0, Ft = 0, Wt = 0, Vt = 0, Jt = 0, Qt = 0, Xt = 0, Yt = 0, $t = 0, te = 0, ee = 0, re = 0, ie = 0, ne = 0, se = 0, ae = 0, he = 0, oe = 0, fe = 0, ue = 0, le = 0, ce = 0, we = 0, ye = 0, pe = 0, ve = 0, de = 0, ge = 0, be = 0, me = 0, Se = 0, Ae = 0, Ee = 0, _e = 0, ke = 0, Ue = 0, xe = 0, Le = 0, He = 0, Ie = 0, qe = 0, Me = 0, Ce = 0, Ze = 0, ze = 0, Oe = 0, Te = 0, Re = 0, Be = 0, Pe = 0, Ke = 0, Ne = 0, je = 0, De = 0, Ge = 0, Fe = 0, We = 0, Ve = 0, Je = 0, Qe = 0, Xe = 0, Ye = 0, $e = 0, tr = 0, er = 0, rr = 0, ir = 0, nr = 0, sr = 0, ar = 0, hr = 0, or = 0, fr = 0, ur = 0, lr = 0, cr = 0, wr = 0, yr = 0, pr = 0, vr = 0, dr = 0, gr = 0, br = 0, mr = 0, Sr = 0, Ar = 0, Er = 0, _r = 0, kr = 0, Ur = 0, xr = 0, Lr = 0, Hr = 0, Ir = 0, qr = 0, Mr = 0, Cr = 0, Zr = 0, zr = 0, Or = 0, Tr = 0, Rr = 0, Br = 0, Pr = 0, Kr = 0, Nr = 0, jr = 0, Dr = 0, Gr = 0, Fr = 0, Wr = 0, Vr = 0, Jr = 0, Qr = 0, Xr = 0, Yr = 0, $r = 0, ti = 0, ei = 0, ri = 0, ii = 0, ni = 0, si = 0, ai = 0, hi = 0, oi = 0, fi = 0, ui = 0, li = 0, ci = 0, wi = 0, yi = 0, pi = 0, vi = 0, di = 0, gi = 0, bi = 0, mi = 0, Si = 0, Ai = 0, Ei = 0, _i = 0, ki = 0, Ui = 0, xi = 0, Li = 0, Hi = 0, Ii = 0, qi = 0, Mi = 0, Ci = 0, Zi = 0, zi = 0, Oi = 0, Ti = 0, Ri = 0, Bi = 0, Pi = 0, Ki = 0, Ni = 0, ji = 0, Di = 0, Gi = 0, Fi = 0, Wi = 0, Vi = 0, Ji = 0, Qi = 0, Xi = 0;
        var Yi = new t.Uint8Array(r);
        function $i() {
            var t = 0;
            D = L ^ Yi[t | K] ^ 1;
            G = H ^ Yi[t | N];
            F = I ^ Yi[t | j];
            W = q ^ Yi[t | P];
            V = M ^ D;
            J = C ^ G;
            Q = Z ^ F;
            X = z ^ W;
            Y = O ^ V;
            $ = T ^ J;
            tt = R ^ Q;
            et = B ^ X;
            rt = P ^ Y;
            it = K ^ $;
            nt = N ^ tt;
            st = j ^ et;
            at = D ^ Yi[t | it] ^ 2;
            ht = G ^ Yi[t | nt];
            ot = F ^ Yi[t | st];
            ft = W ^ Yi[t | rt];
            ut = V ^ at;
            lt = J ^ ht;
            ct = Q ^ ot;
            wt = X ^ ft;
            yt = Y ^ ut;
            pt = $ ^ lt;
            vt = tt ^ ct;
            dt = et ^ wt;
            gt = rt ^ yt;
            bt = it ^ pt;
            mt = nt ^ vt;
            St = st ^ dt;
            At = at ^ Yi[t | bt] ^ 4;
            Et = ht ^ Yi[t | mt];
            _t = ot ^ Yi[t | St];
            kt = ft ^ Yi[t | gt];
            Ut = ut ^ At;
            xt = lt ^ Et;
            Lt = ct ^ _t;
            Ht = wt ^ kt;
            It = yt ^ Ut;
            qt = pt ^ xt;
            Mt = vt ^ Lt;
            Ct = dt ^ Ht;
            Zt = gt ^ It;
            zt = bt ^ qt;
            Ot = mt ^ Mt;
            Tt = St ^ Ct;
            Rt = At ^ Yi[t | zt] ^ 8;
            Bt = Et ^ Yi[t | Ot];
            Pt = _t ^ Yi[t | Tt];
            Kt = kt ^ Yi[t | Zt];
            Nt = Ut ^ Rt;
            jt = xt ^ Bt;
            Dt = Lt ^ Pt;
            Gt = Ht ^ Kt;
            Ft = It ^ Nt;
            Wt = qt ^ jt;
            Vt = Mt ^ Dt;
            Jt = Ct ^ Gt;
            Qt = Zt ^ Ft;
            Xt = zt ^ Wt;
            Yt = Ot ^ Vt;
            $t = Tt ^ Jt;
            te = Rt ^ Yi[t | Xt] ^ 16;
            ee = Bt ^ Yi[t | Yt];
            re = Pt ^ Yi[t | $t];
            ie = Kt ^ Yi[t | Qt];
            ne = Nt ^ te;
            se = jt ^ ee;
            ae = Dt ^ re;
            he = Gt ^ ie;
            oe = Ft ^ ne;
            fe = Wt ^ se;
            ue = Vt ^ ae;
            le = Jt ^ he;
            ce = Qt ^ oe;
            we = Xt ^ fe;
            ye = Yt ^ ue;
            pe = $t ^ le;
            ve = te ^ Yi[t | we] ^ 32;
            de = ee ^ Yi[t | ye];
            ge = re ^ Yi[t | pe];
            be = ie ^ Yi[t | ce];
            me = ne ^ ve;
            Se = se ^ de;
            Ae = ae ^ ge;
            Ee = he ^ be;
            _e = oe ^ me;
            ke = fe ^ Se;
            Ue = ue ^ Ae;
            xe = le ^ Ee;
            Le = ce ^ _e;
            He = we ^ ke;
            Ie = ye ^ Ue;
            qe = pe ^ xe;
            Me = ve ^ Yi[t | He] ^ 64;
            Ce = de ^ Yi[t | Ie];
            Ze = ge ^ Yi[t | qe];
            ze = be ^ Yi[t | Le];
            Oe = me ^ Me;
            Te = Se ^ Ce;
            Re = Ae ^ Ze;
            Be = Ee ^ ze;
            Pe = _e ^ Oe;
            Ke = ke ^ Te;
            Ne = Ue ^ Re;
            je = xe ^ Be;
            De = Le ^ Pe;
            Ge = He ^ Ke;
            Fe = Ie ^ Ne;
            We = qe ^ je;
            Ve = Me ^ Yi[t | Ge] ^ 128;
            Je = Ce ^ Yi[t | Fe];
            Qe = Ze ^ Yi[t | We];
            Xe = ze ^ Yi[t | De];
            Ye = Oe ^ Ve;
            $e = Te ^ Je;
            tr = Re ^ Qe;
            er = Be ^ Xe;
            rr = Pe ^ Ye;
            ir = Ke ^ $e;
            nr = Ne ^ tr;
            sr = je ^ er;
            ar = De ^ rr;
            hr = Ge ^ ir;
            or = Fe ^ nr;
            fr = We ^ sr;
            ur = Ve ^ Yi[t | hr] ^ 27;
            lr = Je ^ Yi[t | or];
            cr = Qe ^ Yi[t | fr];
            wr = Xe ^ Yi[t | ar];
            yr = Ye ^ ur;
            pr = $e ^ lr;
            vr = tr ^ cr;
            dr = er ^ wr;
            gr = rr ^ yr;
            br = ir ^ pr;
            mr = nr ^ vr;
            Sr = sr ^ dr;
            Ar = ar ^ gr;
            Er = hr ^ br;
            _r = or ^ mr;
            kr = fr ^ Sr;
            Ur = ur ^ Yi[t | Er] ^ 54;
            xr = lr ^ Yi[t | _r];
            Lr = cr ^ Yi[t | kr];
            Hr = wr ^ Yi[t | Ar];
            Ir = yr ^ Ur;
            qr = pr ^ xr;
            Mr = vr ^ Lr;
            Cr = dr ^ Hr;
            Zr = gr ^ Ir;
            zr = br ^ qr;
            Or = mr ^ Mr;
            Tr = Sr ^ Cr;
            Rr = Ar ^ Zr;
            Br = Er ^ zr;
            Pr = _r ^ Or;
            Kr = kr ^ Tr;
        }
        function tn() {
            var t = 0;
            at = L ^ Yi[t | it] ^ 1;
            ht = H ^ Yi[t | nt];
            ot = I ^ Yi[t | st];
            ft = q ^ Yi[t | rt];
            ut = M ^ at;
            lt = C ^ ht;
            ct = Z ^ ot;
            wt = z ^ ft;
            yt = O ^ ut;
            pt = T ^ lt;
            vt = R ^ ct;
            dt = B ^ wt;
            gt = P ^ yt;
            bt = K ^ pt;
            mt = N ^ vt;
            St = j ^ dt;
            At = D ^ Yi[t | gt];
            Et = G ^ Yi[t | bt];
            _t = F ^ Yi[t | mt];
            kt = W ^ Yi[t | St];
            Ut = V ^ At;
            xt = J ^ Et;
            Lt = Q ^ _t;
            Ht = X ^ kt;
            It = Y ^ Ut;
            qt = $ ^ xt;
            Mt = tt ^ Lt;
            Ct = et ^ Ht;
            Zt = rt ^ It;
            zt = it ^ qt;
            Ot = nt ^ Mt;
            Tt = st ^ Ct;
            Rt = at ^ Yi[t | zt] ^ 2;
            Bt = ht ^ Yi[t | Ot];
            Pt = ot ^ Yi[t | Tt];
            Kt = ft ^ Yi[t | Zt];
            Nt = ut ^ Rt;
            jt = lt ^ Bt;
            Dt = ct ^ Pt;
            Gt = wt ^ Kt;
            Ft = yt ^ Nt;
            Wt = pt ^ jt;
            Vt = vt ^ Dt;
            Jt = dt ^ Gt;
            Qt = gt ^ Ft;
            Xt = bt ^ Wt;
            Yt = mt ^ Vt;
            $t = St ^ Jt;
            te = At ^ Yi[t | Qt];
            ee = Et ^ Yi[t | Xt];
            re = _t ^ Yi[t | Yt];
            ie = kt ^ Yi[t | $t];
            ne = Ut ^ te;
            se = xt ^ ee;
            ae = Lt ^ re;
            he = Ht ^ ie;
            oe = It ^ ne;
            fe = qt ^ se;
            ue = Mt ^ ae;
            le = Ct ^ he;
            ce = Zt ^ oe;
            we = zt ^ fe;
            ye = Ot ^ ue;
            pe = Tt ^ le;
            ve = Rt ^ Yi[t | we] ^ 4;
            de = Bt ^ Yi[t | ye];
            ge = Pt ^ Yi[t | pe];
            be = Kt ^ Yi[t | ce];
            me = Nt ^ ve;
            Se = jt ^ de;
            Ae = Dt ^ ge;
            Ee = Gt ^ be;
            _e = Ft ^ me;
            ke = Wt ^ Se;
            Ue = Vt ^ Ae;
            xe = Jt ^ Ee;
            Le = Qt ^ _e;
            He = Xt ^ ke;
            Ie = Yt ^ Ue;
            qe = $t ^ xe;
            Me = te ^ Yi[t | Le];
            Ce = ee ^ Yi[t | He];
            Ze = re ^ Yi[t | Ie];
            ze = ie ^ Yi[t | qe];
            Oe = ne ^ Me;
            Te = se ^ Ce;
            Re = ae ^ Ze;
            Be = he ^ ze;
            Pe = oe ^ Oe;
            Ke = fe ^ Te;
            Ne = ue ^ Re;
            je = le ^ Be;
            De = ce ^ Pe;
            Ge = we ^ Ke;
            Fe = ye ^ Ne;
            We = pe ^ je;
            Ve = ve ^ Yi[t | Ge] ^ 8;
            Je = de ^ Yi[t | Fe];
            Qe = ge ^ Yi[t | We];
            Xe = be ^ Yi[t | De];
            Ye = me ^ Ve;
            $e = Se ^ Je;
            tr = Ae ^ Qe;
            er = Ee ^ Xe;
            rr = _e ^ Ye;
            ir = ke ^ $e;
            nr = Ue ^ tr;
            sr = xe ^ er;
            ar = Le ^ rr;
            hr = He ^ ir;
            or = Ie ^ nr;
            fr = qe ^ sr;
            ur = Me ^ Yi[t | ar];
            lr = Ce ^ Yi[t | hr];
            cr = Ze ^ Yi[t | or];
            wr = ze ^ Yi[t | fr];
            yr = Oe ^ ur;
            pr = Te ^ lr;
            vr = Re ^ cr;
            dr = Be ^ wr;
            gr = Pe ^ yr;
            br = Ke ^ pr;
            mr = Ne ^ vr;
            Sr = je ^ dr;
            Ar = De ^ gr;
            Er = Ge ^ br;
            _r = Fe ^ mr;
            kr = We ^ Sr;
            Ur = Ve ^ Yi[t | Er] ^ 16;
            xr = Je ^ Yi[t | _r];
            Lr = Qe ^ Yi[t | kr];
            Hr = Xe ^ Yi[t | Ar];
            Ir = Ye ^ Ur;
            qr = $e ^ xr;
            Mr = tr ^ Lr;
            Cr = er ^ Hr;
            Zr = rr ^ Ir;
            zr = ir ^ qr;
            Or = nr ^ Mr;
            Tr = sr ^ Cr;
            Rr = ar ^ Zr;
            Br = hr ^ zr;
            Pr = or ^ Or;
            Kr = fr ^ Tr;
            Nr = ur ^ Yi[t | Rr];
            jr = lr ^ Yi[t | Br];
            Dr = cr ^ Yi[t | Pr];
            Gr = wr ^ Yi[t | Kr];
            Fr = yr ^ Nr;
            Wr = pr ^ jr;
            Vr = vr ^ Dr;
            Jr = dr ^ Gr;
            Qr = gr ^ Fr;
            Xr = br ^ Wr;
            Yr = mr ^ Vr;
            $r = Sr ^ Jr;
            ti = Ar ^ Qr;
            ei = Er ^ Xr;
            ri = _r ^ Yr;
            ii = kr ^ $r;
            ni = Ur ^ Yi[t | ei] ^ 32;
            si = xr ^ Yi[t | ri];
            ai = Lr ^ Yi[t | ii];
            hi = Hr ^ Yi[t | ti];
            oi = Ir ^ ni;
            fi = qr ^ si;
            ui = Mr ^ ai;
            li = Cr ^ hi;
            ci = Zr ^ oi;
            wi = zr ^ fi;
            yi = Or ^ ui;
            pi = Tr ^ li;
            vi = Rr ^ ci;
            di = Br ^ wi;
            gi = Pr ^ yi;
            bi = Kr ^ pi;
            mi = Nr ^ Yi[t | vi];
            Si = jr ^ Yi[t | di];
            Ai = Dr ^ Yi[t | gi];
            Ei = Gr ^ Yi[t | bi];
            _i = Fr ^ mi;
            ki = Wr ^ Si;
            Ui = Vr ^ Ai;
            xi = Jr ^ Ei;
            Li = Qr ^ _i;
            Hi = Xr ^ ki;
            Ii = Yr ^ Ui;
            qi = $r ^ xi;
            Mi = ti ^ Li;
            Ci = ei ^ Hi;
            Zi = ri ^ Ii;
            zi = ii ^ qi;
            Oi = ni ^ Yi[t | Ci] ^ 64;
            Ti = si ^ Yi[t | Zi];
            Ri = ai ^ Yi[t | zi];
            Bi = hi ^ Yi[t | Mi];
            Pi = oi ^ Oi;
            Ki = fi ^ Ti;
            Ni = ui ^ Ri;
            ji = li ^ Bi;
            Di = ci ^ Pi;
            Gi = wi ^ Ki;
            Fi = yi ^ Ni;
            Wi = pi ^ ji;
            Vi = vi ^ Di;
            Ji = di ^ Gi;
            Qi = gi ^ Fi;
            Xi = bi ^ Wi;
        }
        function en(t, e, r, m, S, A, E, _, k, U, x, $i, tn, en, rn, nn) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            m = m | 0;
            S = S | 0;
            A = A | 0;
            E = E | 0;
            _ = _ | 0;
            k = k | 0;
            U = U | 0;
            x = x | 0;
            $i = $i | 0;
            tn = tn | 0;
            en = en | 0;
            rn = rn | 0;
            nn = nn | 0;
            var sn = 0, an = 0, hn = 0, on = 0, fn = 0, un = 0, ln = 0, cn = 0, wn = 0, yn = 0, pn = 0, vn = 0, dn = 0, gn = 0, bn = 0, mn = 0, Sn = 0, An = 512, En = 768;
            t = t ^ L;
            e = e ^ H;
            r = r ^ I;
            m = m ^ q;
            S = S ^ M;
            A = A ^ C;
            E = E ^ Z;
            _ = _ ^ z;
            k = k ^ O;
            U = U ^ T;
            x = x ^ R;
            $i = $i ^ B;
            tn = tn ^ P;
            en = en ^ K;
            rn = rn ^ N;
            nn = nn ^ j;
            sn = Yi[An | t] ^ Yi[En | A] ^ Yi[Sn | x] ^ Yi[Sn | nn] ^ D;
            an = Yi[Sn | t] ^ Yi[An | A] ^ Yi[En | x] ^ Yi[Sn | nn] ^ G;
            hn = Yi[Sn | t] ^ Yi[Sn | A] ^ Yi[An | x] ^ Yi[En | nn] ^ F;
            on = Yi[En | t] ^ Yi[Sn | A] ^ Yi[Sn | x] ^ Yi[An | nn] ^ W;
            fn = Yi[An | S] ^ Yi[En | U] ^ Yi[Sn | rn] ^ Yi[Sn | m] ^ V;
            un = Yi[Sn | S] ^ Yi[An | U] ^ Yi[En | rn] ^ Yi[Sn | m] ^ J;
            ln = Yi[Sn | S] ^ Yi[Sn | U] ^ Yi[An | rn] ^ Yi[En | m] ^ Q;
            cn = Yi[En | S] ^ Yi[Sn | U] ^ Yi[Sn | rn] ^ Yi[An | m] ^ X;
            wn = Yi[An | k] ^ Yi[En | en] ^ Yi[Sn | r] ^ Yi[Sn | _] ^ Y;
            yn = Yi[Sn | k] ^ Yi[An | en] ^ Yi[En | r] ^ Yi[Sn | _] ^ $;
            pn = Yi[Sn | k] ^ Yi[Sn | en] ^ Yi[An | r] ^ Yi[En | _] ^ tt;
            vn = Yi[En | k] ^ Yi[Sn | en] ^ Yi[Sn | r] ^ Yi[An | _] ^ et;
            dn = Yi[An | tn] ^ Yi[En | e] ^ Yi[Sn | E] ^ Yi[Sn | $i] ^ rt;
            gn = Yi[Sn | tn] ^ Yi[An | e] ^ Yi[En | E] ^ Yi[Sn | $i] ^ it;
            bn = Yi[Sn | tn] ^ Yi[Sn | e] ^ Yi[An | E] ^ Yi[En | $i] ^ nt;
            mn = Yi[En | tn] ^ Yi[Sn | e] ^ Yi[Sn | E] ^ Yi[An | $i] ^ st;
            t = Yi[An | sn] ^ Yi[En | un] ^ Yi[Sn | pn] ^ Yi[Sn | mn] ^ at;
            e = Yi[Sn | sn] ^ Yi[An | un] ^ Yi[En | pn] ^ Yi[Sn | mn] ^ ht;
            r = Yi[Sn | sn] ^ Yi[Sn | un] ^ Yi[An | pn] ^ Yi[En | mn] ^ ot;
            m = Yi[En | sn] ^ Yi[Sn | un] ^ Yi[Sn | pn] ^ Yi[An | mn] ^ ft;
            S = Yi[An | fn] ^ Yi[En | yn] ^ Yi[Sn | bn] ^ Yi[Sn | on] ^ ut;
            A = Yi[Sn | fn] ^ Yi[An | yn] ^ Yi[En | bn] ^ Yi[Sn | on] ^ lt;
            E = Yi[Sn | fn] ^ Yi[Sn | yn] ^ Yi[An | bn] ^ Yi[En | on] ^ ct;
            _ = Yi[En | fn] ^ Yi[Sn | yn] ^ Yi[Sn | bn] ^ Yi[An | on] ^ wt;
            k = Yi[An | wn] ^ Yi[En | gn] ^ Yi[Sn | hn] ^ Yi[Sn | cn] ^ yt;
            U = Yi[Sn | wn] ^ Yi[An | gn] ^ Yi[En | hn] ^ Yi[Sn | cn] ^ pt;
            x = Yi[Sn | wn] ^ Yi[Sn | gn] ^ Yi[An | hn] ^ Yi[En | cn] ^ vt;
            $i = Yi[En | wn] ^ Yi[Sn | gn] ^ Yi[Sn | hn] ^ Yi[An | cn] ^ dt;
            tn = Yi[An | dn] ^ Yi[En | an] ^ Yi[Sn | ln] ^ Yi[Sn | vn] ^ gt;
            en = Yi[Sn | dn] ^ Yi[An | an] ^ Yi[En | ln] ^ Yi[Sn | vn] ^ bt;
            rn = Yi[Sn | dn] ^ Yi[Sn | an] ^ Yi[An | ln] ^ Yi[En | vn] ^ mt;
            nn = Yi[En | dn] ^ Yi[Sn | an] ^ Yi[Sn | ln] ^ Yi[An | vn] ^ St;
            sn = Yi[An | t] ^ Yi[En | A] ^ Yi[Sn | x] ^ Yi[Sn | nn] ^ At;
            an = Yi[Sn | t] ^ Yi[An | A] ^ Yi[En | x] ^ Yi[Sn | nn] ^ Et;
            hn = Yi[Sn | t] ^ Yi[Sn | A] ^ Yi[An | x] ^ Yi[En | nn] ^ _t;
            on = Yi[En | t] ^ Yi[Sn | A] ^ Yi[Sn | x] ^ Yi[An | nn] ^ kt;
            fn = Yi[An | S] ^ Yi[En | U] ^ Yi[Sn | rn] ^ Yi[Sn | m] ^ Ut;
            un = Yi[Sn | S] ^ Yi[An | U] ^ Yi[En | rn] ^ Yi[Sn | m] ^ xt;
            ln = Yi[Sn | S] ^ Yi[Sn | U] ^ Yi[An | rn] ^ Yi[En | m] ^ Lt;
            cn = Yi[En | S] ^ Yi[Sn | U] ^ Yi[Sn | rn] ^ Yi[An | m] ^ Ht;
            wn = Yi[An | k] ^ Yi[En | en] ^ Yi[Sn | r] ^ Yi[Sn | _] ^ It;
            yn = Yi[Sn | k] ^ Yi[An | en] ^ Yi[En | r] ^ Yi[Sn | _] ^ qt;
            pn = Yi[Sn | k] ^ Yi[Sn | en] ^ Yi[An | r] ^ Yi[En | _] ^ Mt;
            vn = Yi[En | k] ^ Yi[Sn | en] ^ Yi[Sn | r] ^ Yi[An | _] ^ Ct;
            dn = Yi[An | tn] ^ Yi[En | e] ^ Yi[Sn | E] ^ Yi[Sn | $i] ^ Zt;
            gn = Yi[Sn | tn] ^ Yi[An | e] ^ Yi[En | E] ^ Yi[Sn | $i] ^ zt;
            bn = Yi[Sn | tn] ^ Yi[Sn | e] ^ Yi[An | E] ^ Yi[En | $i] ^ Ot;
            mn = Yi[En | tn] ^ Yi[Sn | e] ^ Yi[Sn | E] ^ Yi[An | $i] ^ Tt;
            t = Yi[An | sn] ^ Yi[En | un] ^ Yi[Sn | pn] ^ Yi[Sn | mn] ^ Rt;
            e = Yi[Sn | sn] ^ Yi[An | un] ^ Yi[En | pn] ^ Yi[Sn | mn] ^ Bt;
            r = Yi[Sn | sn] ^ Yi[Sn | un] ^ Yi[An | pn] ^ Yi[En | mn] ^ Pt;
            m = Yi[En | sn] ^ Yi[Sn | un] ^ Yi[Sn | pn] ^ Yi[An | mn] ^ Kt;
            S = Yi[An | fn] ^ Yi[En | yn] ^ Yi[Sn | bn] ^ Yi[Sn | on] ^ Nt;
            A = Yi[Sn | fn] ^ Yi[An | yn] ^ Yi[En | bn] ^ Yi[Sn | on] ^ jt;
            E = Yi[Sn | fn] ^ Yi[Sn | yn] ^ Yi[An | bn] ^ Yi[En | on] ^ Dt;
            _ = Yi[En | fn] ^ Yi[Sn | yn] ^ Yi[Sn | bn] ^ Yi[An | on] ^ Gt;
            k = Yi[An | wn] ^ Yi[En | gn] ^ Yi[Sn | hn] ^ Yi[Sn | cn] ^ Ft;
            U = Yi[Sn | wn] ^ Yi[An | gn] ^ Yi[En | hn] ^ Yi[Sn | cn] ^ Wt;
            x = Yi[Sn | wn] ^ Yi[Sn | gn] ^ Yi[An | hn] ^ Yi[En | cn] ^ Vt;
            $i = Yi[En | wn] ^ Yi[Sn | gn] ^ Yi[Sn | hn] ^ Yi[An | cn] ^ Jt;
            tn = Yi[An | dn] ^ Yi[En | an] ^ Yi[Sn | ln] ^ Yi[Sn | vn] ^ Qt;
            en = Yi[Sn | dn] ^ Yi[An | an] ^ Yi[En | ln] ^ Yi[Sn | vn] ^ Xt;
            rn = Yi[Sn | dn] ^ Yi[Sn | an] ^ Yi[An | ln] ^ Yi[En | vn] ^ Yt;
            nn = Yi[En | dn] ^ Yi[Sn | an] ^ Yi[Sn | ln] ^ Yi[An | vn] ^ $t;
            sn = Yi[An | t] ^ Yi[En | A] ^ Yi[Sn | x] ^ Yi[Sn | nn] ^ te;
            an = Yi[Sn | t] ^ Yi[An | A] ^ Yi[En | x] ^ Yi[Sn | nn] ^ ee;
            hn = Yi[Sn | t] ^ Yi[Sn | A] ^ Yi[An | x] ^ Yi[En | nn] ^ re;
            on = Yi[En | t] ^ Yi[Sn | A] ^ Yi[Sn | x] ^ Yi[An | nn] ^ ie;
            fn = Yi[An | S] ^ Yi[En | U] ^ Yi[Sn | rn] ^ Yi[Sn | m] ^ ne;
            un = Yi[Sn | S] ^ Yi[An | U] ^ Yi[En | rn] ^ Yi[Sn | m] ^ se;
            ln = Yi[Sn | S] ^ Yi[Sn | U] ^ Yi[An | rn] ^ Yi[En | m] ^ ae;
            cn = Yi[En | S] ^ Yi[Sn | U] ^ Yi[Sn | rn] ^ Yi[An | m] ^ he;
            wn = Yi[An | k] ^ Yi[En | en] ^ Yi[Sn | r] ^ Yi[Sn | _] ^ oe;
            yn = Yi[Sn | k] ^ Yi[An | en] ^ Yi[En | r] ^ Yi[Sn | _] ^ fe;
            pn = Yi[Sn | k] ^ Yi[Sn | en] ^ Yi[An | r] ^ Yi[En | _] ^ ue;
            vn = Yi[En | k] ^ Yi[Sn | en] ^ Yi[Sn | r] ^ Yi[An | _] ^ le;
            dn = Yi[An | tn] ^ Yi[En | e] ^ Yi[Sn | E] ^ Yi[Sn | $i] ^ ce;
            gn = Yi[Sn | tn] ^ Yi[An | e] ^ Yi[En | E] ^ Yi[Sn | $i] ^ we;
            bn = Yi[Sn | tn] ^ Yi[Sn | e] ^ Yi[An | E] ^ Yi[En | $i] ^ ye;
            mn = Yi[En | tn] ^ Yi[Sn | e] ^ Yi[Sn | E] ^ Yi[An | $i] ^ pe;
            t = Yi[An | sn] ^ Yi[En | un] ^ Yi[Sn | pn] ^ Yi[Sn | mn] ^ ve;
            e = Yi[Sn | sn] ^ Yi[An | un] ^ Yi[En | pn] ^ Yi[Sn | mn] ^ de;
            r = Yi[Sn | sn] ^ Yi[Sn | un] ^ Yi[An | pn] ^ Yi[En | mn] ^ ge;
            m = Yi[En | sn] ^ Yi[Sn | un] ^ Yi[Sn | pn] ^ Yi[An | mn] ^ be;
            S = Yi[An | fn] ^ Yi[En | yn] ^ Yi[Sn | bn] ^ Yi[Sn | on] ^ me;
            A = Yi[Sn | fn] ^ Yi[An | yn] ^ Yi[En | bn] ^ Yi[Sn | on] ^ Se;
            E = Yi[Sn | fn] ^ Yi[Sn | yn] ^ Yi[An | bn] ^ Yi[En | on] ^ Ae;
            _ = Yi[En | fn] ^ Yi[Sn | yn] ^ Yi[Sn | bn] ^ Yi[An | on] ^ Ee;
            k = Yi[An | wn] ^ Yi[En | gn] ^ Yi[Sn | hn] ^ Yi[Sn | cn] ^ _e;
            U = Yi[Sn | wn] ^ Yi[An | gn] ^ Yi[En | hn] ^ Yi[Sn | cn] ^ ke;
            x = Yi[Sn | wn] ^ Yi[Sn | gn] ^ Yi[An | hn] ^ Yi[En | cn] ^ Ue;
            $i = Yi[En | wn] ^ Yi[Sn | gn] ^ Yi[Sn | hn] ^ Yi[An | cn] ^ xe;
            tn = Yi[An | dn] ^ Yi[En | an] ^ Yi[Sn | ln] ^ Yi[Sn | vn] ^ Le;
            en = Yi[Sn | dn] ^ Yi[An | an] ^ Yi[En | ln] ^ Yi[Sn | vn] ^ He;
            rn = Yi[Sn | dn] ^ Yi[Sn | an] ^ Yi[An | ln] ^ Yi[En | vn] ^ Ie;
            nn = Yi[En | dn] ^ Yi[Sn | an] ^ Yi[Sn | ln] ^ Yi[An | vn] ^ qe;
            sn = Yi[An | t] ^ Yi[En | A] ^ Yi[Sn | x] ^ Yi[Sn | nn] ^ Me;
            an = Yi[Sn | t] ^ Yi[An | A] ^ Yi[En | x] ^ Yi[Sn | nn] ^ Ce;
            hn = Yi[Sn | t] ^ Yi[Sn | A] ^ Yi[An | x] ^ Yi[En | nn] ^ Ze;
            on = Yi[En | t] ^ Yi[Sn | A] ^ Yi[Sn | x] ^ Yi[An | nn] ^ ze;
            fn = Yi[An | S] ^ Yi[En | U] ^ Yi[Sn | rn] ^ Yi[Sn | m] ^ Oe;
            un = Yi[Sn | S] ^ Yi[An | U] ^ Yi[En | rn] ^ Yi[Sn | m] ^ Te;
            ln = Yi[Sn | S] ^ Yi[Sn | U] ^ Yi[An | rn] ^ Yi[En | m] ^ Re;
            cn = Yi[En | S] ^ Yi[Sn | U] ^ Yi[Sn | rn] ^ Yi[An | m] ^ Be;
            wn = Yi[An | k] ^ Yi[En | en] ^ Yi[Sn | r] ^ Yi[Sn | _] ^ Pe;
            yn = Yi[Sn | k] ^ Yi[An | en] ^ Yi[En | r] ^ Yi[Sn | _] ^ Ke;
            pn = Yi[Sn | k] ^ Yi[Sn | en] ^ Yi[An | r] ^ Yi[En | _] ^ Ne;
            vn = Yi[En | k] ^ Yi[Sn | en] ^ Yi[Sn | r] ^ Yi[An | _] ^ je;
            dn = Yi[An | tn] ^ Yi[En | e] ^ Yi[Sn | E] ^ Yi[Sn | $i] ^ De;
            gn = Yi[Sn | tn] ^ Yi[An | e] ^ Yi[En | E] ^ Yi[Sn | $i] ^ Ge;
            bn = Yi[Sn | tn] ^ Yi[Sn | e] ^ Yi[An | E] ^ Yi[En | $i] ^ Fe;
            mn = Yi[En | tn] ^ Yi[Sn | e] ^ Yi[Sn | E] ^ Yi[An | $i] ^ We;
            t = Yi[An | sn] ^ Yi[En | un] ^ Yi[Sn | pn] ^ Yi[Sn | mn] ^ Ve;
            e = Yi[Sn | sn] ^ Yi[An | un] ^ Yi[En | pn] ^ Yi[Sn | mn] ^ Je;
            r = Yi[Sn | sn] ^ Yi[Sn | un] ^ Yi[An | pn] ^ Yi[En | mn] ^ Qe;
            m = Yi[En | sn] ^ Yi[Sn | un] ^ Yi[Sn | pn] ^ Yi[An | mn] ^ Xe;
            S = Yi[An | fn] ^ Yi[En | yn] ^ Yi[Sn | bn] ^ Yi[Sn | on] ^ Ye;
            A = Yi[Sn | fn] ^ Yi[An | yn] ^ Yi[En | bn] ^ Yi[Sn | on] ^ $e;
            E = Yi[Sn | fn] ^ Yi[Sn | yn] ^ Yi[An | bn] ^ Yi[En | on] ^ tr;
            _ = Yi[En | fn] ^ Yi[Sn | yn] ^ Yi[Sn | bn] ^ Yi[An | on] ^ er;
            k = Yi[An | wn] ^ Yi[En | gn] ^ Yi[Sn | hn] ^ Yi[Sn | cn] ^ rr;
            U = Yi[Sn | wn] ^ Yi[An | gn] ^ Yi[En | hn] ^ Yi[Sn | cn] ^ ir;
            x = Yi[Sn | wn] ^ Yi[Sn | gn] ^ Yi[An | hn] ^ Yi[En | cn] ^ nr;
            $i = Yi[En | wn] ^ Yi[Sn | gn] ^ Yi[Sn | hn] ^ Yi[An | cn] ^ sr;
            tn = Yi[An | dn] ^ Yi[En | an] ^ Yi[Sn | ln] ^ Yi[Sn | vn] ^ ar;
            en = Yi[Sn | dn] ^ Yi[An | an] ^ Yi[En | ln] ^ Yi[Sn | vn] ^ hr;
            rn = Yi[Sn | dn] ^ Yi[Sn | an] ^ Yi[An | ln] ^ Yi[En | vn] ^ or;
            nn = Yi[En | dn] ^ Yi[Sn | an] ^ Yi[Sn | ln] ^ Yi[An | vn] ^ fr;
            sn = Yi[An | t] ^ Yi[En | A] ^ Yi[Sn | x] ^ Yi[Sn | nn] ^ ur;
            an = Yi[Sn | t] ^ Yi[An | A] ^ Yi[En | x] ^ Yi[Sn | nn] ^ lr;
            hn = Yi[Sn | t] ^ Yi[Sn | A] ^ Yi[An | x] ^ Yi[En | nn] ^ cr;
            on = Yi[En | t] ^ Yi[Sn | A] ^ Yi[Sn | x] ^ Yi[An | nn] ^ wr;
            fn = Yi[An | S] ^ Yi[En | U] ^ Yi[Sn | rn] ^ Yi[Sn | m] ^ yr;
            un = Yi[Sn | S] ^ Yi[An | U] ^ Yi[En | rn] ^ Yi[Sn | m] ^ pr;
            ln = Yi[Sn | S] ^ Yi[Sn | U] ^ Yi[An | rn] ^ Yi[En | m] ^ vr;
            cn = Yi[En | S] ^ Yi[Sn | U] ^ Yi[Sn | rn] ^ Yi[An | m] ^ dr;
            wn = Yi[An | k] ^ Yi[En | en] ^ Yi[Sn | r] ^ Yi[Sn | _] ^ gr;
            yn = Yi[Sn | k] ^ Yi[An | en] ^ Yi[En | r] ^ Yi[Sn | _] ^ br;
            pn = Yi[Sn | k] ^ Yi[Sn | en] ^ Yi[An | r] ^ Yi[En | _] ^ mr;
            vn = Yi[En | k] ^ Yi[Sn | en] ^ Yi[Sn | r] ^ Yi[An | _] ^ Sr;
            dn = Yi[An | tn] ^ Yi[En | e] ^ Yi[Sn | E] ^ Yi[Sn | $i] ^ Ar;
            gn = Yi[Sn | tn] ^ Yi[An | e] ^ Yi[En | E] ^ Yi[Sn | $i] ^ Er;
            bn = Yi[Sn | tn] ^ Yi[Sn | e] ^ Yi[An | E] ^ Yi[En | $i] ^ _r;
            mn = Yi[En | tn] ^ Yi[Sn | e] ^ Yi[Sn | E] ^ Yi[An | $i] ^ kr;
            if ((b | 0) == 16) {
                i = Yi[Sn | sn] ^ Ur;
                n = Yi[Sn | un] ^ xr;
                s = Yi[Sn | pn] ^ Lr;
                a = Yi[Sn | mn] ^ Hr;
                h = Yi[Sn | fn] ^ Ir;
                o = Yi[Sn | yn] ^ qr;
                f = Yi[Sn | bn] ^ Mr;
                u = Yi[Sn | on] ^ Cr;
                l = Yi[Sn | wn] ^ Zr;
                c = Yi[Sn | gn] ^ zr;
                w = Yi[Sn | hn] ^ Or;
                y = Yi[Sn | cn] ^ Tr;
                p = Yi[Sn | dn] ^ Rr;
                v = Yi[Sn | an] ^ Br;
                d = Yi[Sn | ln] ^ Pr;
                g = Yi[Sn | vn] ^ Kr;
                return;
            }
            t = Yi[An | sn] ^ Yi[En | un] ^ Yi[Sn | pn] ^ Yi[Sn | mn] ^ Ur;
            e = Yi[Sn | sn] ^ Yi[An | un] ^ Yi[En | pn] ^ Yi[Sn | mn] ^ xr;
            r = Yi[Sn | sn] ^ Yi[Sn | un] ^ Yi[An | pn] ^ Yi[En | mn] ^ Lr;
            m = Yi[En | sn] ^ Yi[Sn | un] ^ Yi[Sn | pn] ^ Yi[An | mn] ^ Hr;
            S = Yi[An | fn] ^ Yi[En | yn] ^ Yi[Sn | bn] ^ Yi[Sn | on] ^ Ir;
            A = Yi[Sn | fn] ^ Yi[An | yn] ^ Yi[En | bn] ^ Yi[Sn | on] ^ qr;
            E = Yi[Sn | fn] ^ Yi[Sn | yn] ^ Yi[An | bn] ^ Yi[En | on] ^ Mr;
            _ = Yi[En | fn] ^ Yi[Sn | yn] ^ Yi[Sn | bn] ^ Yi[An | on] ^ Cr;
            k = Yi[An | wn] ^ Yi[En | gn] ^ Yi[Sn | hn] ^ Yi[Sn | cn] ^ Zr;
            U = Yi[Sn | wn] ^ Yi[An | gn] ^ Yi[En | hn] ^ Yi[Sn | cn] ^ zr;
            x = Yi[Sn | wn] ^ Yi[Sn | gn] ^ Yi[An | hn] ^ Yi[En | cn] ^ Or;
            $i = Yi[En | wn] ^ Yi[Sn | gn] ^ Yi[Sn | hn] ^ Yi[An | cn] ^ Tr;
            tn = Yi[An | dn] ^ Yi[En | an] ^ Yi[Sn | ln] ^ Yi[Sn | vn] ^ Rr;
            en = Yi[Sn | dn] ^ Yi[An | an] ^ Yi[En | ln] ^ Yi[Sn | vn] ^ Br;
            rn = Yi[Sn | dn] ^ Yi[Sn | an] ^ Yi[An | ln] ^ Yi[En | vn] ^ Pr;
            nn = Yi[En | dn] ^ Yi[Sn | an] ^ Yi[Sn | ln] ^ Yi[An | vn] ^ Kr;
            sn = Yi[An | t] ^ Yi[En | A] ^ Yi[Sn | x] ^ Yi[Sn | nn] ^ Nr;
            an = Yi[Sn | t] ^ Yi[An | A] ^ Yi[En | x] ^ Yi[Sn | nn] ^ jr;
            hn = Yi[Sn | t] ^ Yi[Sn | A] ^ Yi[An | x] ^ Yi[En | nn] ^ Dr;
            on = Yi[En | t] ^ Yi[Sn | A] ^ Yi[Sn | x] ^ Yi[An | nn] ^ Gr;
            fn = Yi[An | S] ^ Yi[En | U] ^ Yi[Sn | rn] ^ Yi[Sn | m] ^ Fr;
            un = Yi[Sn | S] ^ Yi[An | U] ^ Yi[En | rn] ^ Yi[Sn | m] ^ Wr;
            ln = Yi[Sn | S] ^ Yi[Sn | U] ^ Yi[An | rn] ^ Yi[En | m] ^ Vr;
            cn = Yi[En | S] ^ Yi[Sn | U] ^ Yi[Sn | rn] ^ Yi[An | m] ^ Jr;
            wn = Yi[An | k] ^ Yi[En | en] ^ Yi[Sn | r] ^ Yi[Sn | _] ^ Qr;
            yn = Yi[Sn | k] ^ Yi[An | en] ^ Yi[En | r] ^ Yi[Sn | _] ^ Xr;
            pn = Yi[Sn | k] ^ Yi[Sn | en] ^ Yi[An | r] ^ Yi[En | _] ^ Yr;
            vn = Yi[En | k] ^ Yi[Sn | en] ^ Yi[Sn | r] ^ Yi[An | _] ^ $r;
            dn = Yi[An | tn] ^ Yi[En | e] ^ Yi[Sn | E] ^ Yi[Sn | $i] ^ ti;
            gn = Yi[Sn | tn] ^ Yi[An | e] ^ Yi[En | E] ^ Yi[Sn | $i] ^ ei;
            bn = Yi[Sn | tn] ^ Yi[Sn | e] ^ Yi[An | E] ^ Yi[En | $i] ^ ri;
            mn = Yi[En | tn] ^ Yi[Sn | e] ^ Yi[Sn | E] ^ Yi[An | $i] ^ ii;
            t = Yi[An | sn] ^ Yi[En | un] ^ Yi[Sn | pn] ^ Yi[Sn | mn] ^ ni;
            e = Yi[Sn | sn] ^ Yi[An | un] ^ Yi[En | pn] ^ Yi[Sn | mn] ^ si;
            r = Yi[Sn | sn] ^ Yi[Sn | un] ^ Yi[An | pn] ^ Yi[En | mn] ^ ai;
            m = Yi[En | sn] ^ Yi[Sn | un] ^ Yi[Sn | pn] ^ Yi[An | mn] ^ hi;
            S = Yi[An | fn] ^ Yi[En | yn] ^ Yi[Sn | bn] ^ Yi[Sn | on] ^ oi;
            A = Yi[Sn | fn] ^ Yi[An | yn] ^ Yi[En | bn] ^ Yi[Sn | on] ^ fi;
            E = Yi[Sn | fn] ^ Yi[Sn | yn] ^ Yi[An | bn] ^ Yi[En | on] ^ ui;
            _ = Yi[En | fn] ^ Yi[Sn | yn] ^ Yi[Sn | bn] ^ Yi[An | on] ^ li;
            k = Yi[An | wn] ^ Yi[En | gn] ^ Yi[Sn | hn] ^ Yi[Sn | cn] ^ ci;
            U = Yi[Sn | wn] ^ Yi[An | gn] ^ Yi[En | hn] ^ Yi[Sn | cn] ^ wi;
            x = Yi[Sn | wn] ^ Yi[Sn | gn] ^ Yi[An | hn] ^ Yi[En | cn] ^ yi;
            $i = Yi[En | wn] ^ Yi[Sn | gn] ^ Yi[Sn | hn] ^ Yi[An | cn] ^ pi;
            tn = Yi[An | dn] ^ Yi[En | an] ^ Yi[Sn | ln] ^ Yi[Sn | vn] ^ vi;
            en = Yi[Sn | dn] ^ Yi[An | an] ^ Yi[En | ln] ^ Yi[Sn | vn] ^ di;
            rn = Yi[Sn | dn] ^ Yi[Sn | an] ^ Yi[An | ln] ^ Yi[En | vn] ^ gi;
            nn = Yi[En | dn] ^ Yi[Sn | an] ^ Yi[Sn | ln] ^ Yi[An | vn] ^ bi;
            sn = Yi[An | t] ^ Yi[En | A] ^ Yi[Sn | x] ^ Yi[Sn | nn] ^ mi;
            an = Yi[Sn | t] ^ Yi[An | A] ^ Yi[En | x] ^ Yi[Sn | nn] ^ Si;
            hn = Yi[Sn | t] ^ Yi[Sn | A] ^ Yi[An | x] ^ Yi[En | nn] ^ Ai;
            on = Yi[En | t] ^ Yi[Sn | A] ^ Yi[Sn | x] ^ Yi[An | nn] ^ Ei;
            fn = Yi[An | S] ^ Yi[En | U] ^ Yi[Sn | rn] ^ Yi[Sn | m] ^ _i;
            un = Yi[Sn | S] ^ Yi[An | U] ^ Yi[En | rn] ^ Yi[Sn | m] ^ ki;
            ln = Yi[Sn | S] ^ Yi[Sn | U] ^ Yi[An | rn] ^ Yi[En | m] ^ Ui;
            cn = Yi[En | S] ^ Yi[Sn | U] ^ Yi[Sn | rn] ^ Yi[An | m] ^ xi;
            wn = Yi[An | k] ^ Yi[En | en] ^ Yi[Sn | r] ^ Yi[Sn | _] ^ Li;
            yn = Yi[Sn | k] ^ Yi[An | en] ^ Yi[En | r] ^ Yi[Sn | _] ^ Hi;
            pn = Yi[Sn | k] ^ Yi[Sn | en] ^ Yi[An | r] ^ Yi[En | _] ^ Ii;
            vn = Yi[En | k] ^ Yi[Sn | en] ^ Yi[Sn | r] ^ Yi[An | _] ^ qi;
            dn = Yi[An | tn] ^ Yi[En | e] ^ Yi[Sn | E] ^ Yi[Sn | $i] ^ Mi;
            gn = Yi[Sn | tn] ^ Yi[An | e] ^ Yi[En | E] ^ Yi[Sn | $i] ^ Ci;
            bn = Yi[Sn | tn] ^ Yi[Sn | e] ^ Yi[An | E] ^ Yi[En | $i] ^ Zi;
            mn = Yi[En | tn] ^ Yi[Sn | e] ^ Yi[Sn | E] ^ Yi[An | $i] ^ zi;
            i = Yi[Sn | sn] ^ Oi;
            n = Yi[Sn | un] ^ Ti;
            s = Yi[Sn | pn] ^ Ri;
            a = Yi[Sn | mn] ^ Bi;
            h = Yi[Sn | fn] ^ Pi;
            o = Yi[Sn | yn] ^ Ki;
            f = Yi[Sn | bn] ^ Ni;
            u = Yi[Sn | on] ^ ji;
            l = Yi[Sn | wn] ^ Di;
            c = Yi[Sn | gn] ^ Gi;
            w = Yi[Sn | hn] ^ Fi;
            y = Yi[Sn | cn] ^ Wi;
            p = Yi[Sn | dn] ^ Vi;
            v = Yi[Sn | an] ^ Ji;
            d = Yi[Sn | ln] ^ Qi;
            g = Yi[Sn | vn] ^ Xi;
        }
        function rn(t, e, r, m, S, A, E, _, k, U, x, $i, tn, en, rn, nn) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            m = m | 0;
            S = S | 0;
            A = A | 0;
            E = E | 0;
            _ = _ | 0;
            k = k | 0;
            U = U | 0;
            x = x | 0;
            $i = $i | 0;
            tn = tn | 0;
            en = en | 0;
            rn = rn | 0;
            nn = nn | 0;
            var sn = 0, an = 0, hn = 0, on = 0, fn = 0, un = 0, ln = 0, cn = 0, wn = 0, yn = 0, pn = 0, vn = 0, dn = 0, gn = 0, bn = 0, mn = 0, Sn = 256, An = 1024, En = 1280, _n = 1536, kn = 1792;
            if ((b | 0) == 32) {
                sn = Yi[Sn | t ^ Oi] ^ mi;
                an = Yi[Sn | en ^ Ji] ^ Si;
                hn = Yi[Sn | x ^ Fi] ^ Ai;
                on = Yi[Sn | _ ^ ji] ^ Ei;
                fn = Yi[Sn | S ^ Pi] ^ _i;
                un = Yi[Sn | e ^ Ti] ^ ki;
                ln = Yi[Sn | rn ^ Qi] ^ Ui;
                cn = Yi[Sn | $i ^ Wi] ^ xi;
                wn = Yi[Sn | k ^ Di] ^ Li;
                yn = Yi[Sn | A ^ Ki] ^ Hi;
                pn = Yi[Sn | r ^ Ri] ^ Ii;
                vn = Yi[Sn | nn ^ Xi] ^ qi;
                dn = Yi[Sn | tn ^ Vi] ^ Mi;
                gn = Yi[Sn | U ^ Gi] ^ Ci;
                bn = Yi[Sn | E ^ Ni] ^ Zi;
                mn = Yi[Sn | m ^ Bi] ^ zi;
                t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
                e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
                r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
                m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
                S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
                A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
                E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
                _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
                k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
                U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
                x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
                $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
                tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
                en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
                rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
                nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
                sn = Yi[Sn | t] ^ ni;
                an = Yi[Sn | e] ^ si;
                hn = Yi[Sn | r] ^ ai;
                on = Yi[Sn | m] ^ hi;
                fn = Yi[Sn | S] ^ oi;
                un = Yi[Sn | A] ^ fi;
                ln = Yi[Sn | E] ^ ui;
                cn = Yi[Sn | _] ^ li;
                wn = Yi[Sn | k] ^ ci;
                yn = Yi[Sn | U] ^ wi;
                pn = Yi[Sn | x] ^ yi;
                vn = Yi[Sn | $i] ^ pi;
                dn = Yi[Sn | tn] ^ vi;
                gn = Yi[Sn | en] ^ di;
                bn = Yi[Sn | rn] ^ gi;
                mn = Yi[Sn | nn] ^ bi;
                t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
                e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
                r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
                m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
                S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
                A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
                E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
                _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
                k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
                U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
                x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
                $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
                tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
                en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
                rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
                nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
                sn = Yi[Sn | t] ^ Nr;
                an = Yi[Sn | e] ^ jr;
                hn = Yi[Sn | r] ^ Dr;
                on = Yi[Sn | m] ^ Gr;
                fn = Yi[Sn | S] ^ Fr;
                un = Yi[Sn | A] ^ Wr;
                ln = Yi[Sn | E] ^ Vr;
                cn = Yi[Sn | _] ^ Jr;
                wn = Yi[Sn | k] ^ Qr;
                yn = Yi[Sn | U] ^ Xr;
                pn = Yi[Sn | x] ^ Yr;
                vn = Yi[Sn | $i] ^ $r;
                dn = Yi[Sn | tn] ^ ti;
                gn = Yi[Sn | en] ^ ei;
                bn = Yi[Sn | rn] ^ ri;
                mn = Yi[Sn | nn] ^ ii;
                t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
                e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
                r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
                m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
                S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
                A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
                E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
                _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
                k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
                U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
                x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
                $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
                tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
                en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
                rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
                nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
                sn = Yi[Sn | t] ^ Ur;
                an = Yi[Sn | e] ^ xr;
                hn = Yi[Sn | r] ^ Lr;
                on = Yi[Sn | m] ^ Hr;
                fn = Yi[Sn | S] ^ Ir;
                un = Yi[Sn | A] ^ qr;
                ln = Yi[Sn | E] ^ Mr;
                cn = Yi[Sn | _] ^ Cr;
                wn = Yi[Sn | k] ^ Zr;
                yn = Yi[Sn | U] ^ zr;
                pn = Yi[Sn | x] ^ Or;
                vn = Yi[Sn | $i] ^ Tr;
                dn = Yi[Sn | tn] ^ Rr;
                gn = Yi[Sn | en] ^ Br;
                bn = Yi[Sn | rn] ^ Pr;
                mn = Yi[Sn | nn] ^ Kr;
                t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
                e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
                r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
                m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
                S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
                A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
                E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
                _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
                k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
                U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
                x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
                $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
                tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
                en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
                rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
                nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
                sn = Yi[Sn | t] ^ ur;
                an = Yi[Sn | e] ^ lr;
                hn = Yi[Sn | r] ^ cr;
                on = Yi[Sn | m] ^ wr;
                fn = Yi[Sn | S] ^ yr;
                un = Yi[Sn | A] ^ pr;
                ln = Yi[Sn | E] ^ vr;
                cn = Yi[Sn | _] ^ dr;
                wn = Yi[Sn | k] ^ gr;
                yn = Yi[Sn | U] ^ br;
                pn = Yi[Sn | x] ^ mr;
                vn = Yi[Sn | $i] ^ Sr;
                dn = Yi[Sn | tn] ^ Ar;
                gn = Yi[Sn | en] ^ Er;
                bn = Yi[Sn | rn] ^ _r;
                mn = Yi[Sn | nn] ^ kr;
            } else {
                sn = Yi[Sn | t ^ Ur] ^ ur;
                an = Yi[Sn | en ^ Br] ^ lr;
                hn = Yi[Sn | x ^ Or] ^ cr;
                on = Yi[Sn | _ ^ Cr] ^ wr;
                fn = Yi[Sn | S ^ Ir] ^ yr;
                un = Yi[Sn | e ^ xr] ^ pr;
                ln = Yi[Sn | rn ^ Pr] ^ vr;
                cn = Yi[Sn | $i ^ Tr] ^ dr;
                wn = Yi[Sn | k ^ Zr] ^ gr;
                yn = Yi[Sn | A ^ qr] ^ br;
                pn = Yi[Sn | r ^ Lr] ^ mr;
                vn = Yi[Sn | nn ^ Kr] ^ Sr;
                dn = Yi[Sn | tn ^ Rr] ^ Ar;
                gn = Yi[Sn | U ^ zr] ^ Er;
                bn = Yi[Sn | E ^ Mr] ^ _r;
                mn = Yi[Sn | m ^ Hr] ^ kr;
            }
            t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
            e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
            r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
            m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
            S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
            A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
            E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
            _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
            k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
            U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
            x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
            $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
            tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
            en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
            rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
            nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
            sn = Yi[Sn | t] ^ Ve;
            an = Yi[Sn | e] ^ Je;
            hn = Yi[Sn | r] ^ Qe;
            on = Yi[Sn | m] ^ Xe;
            fn = Yi[Sn | S] ^ Ye;
            un = Yi[Sn | A] ^ $e;
            ln = Yi[Sn | E] ^ tr;
            cn = Yi[Sn | _] ^ er;
            wn = Yi[Sn | k] ^ rr;
            yn = Yi[Sn | U] ^ ir;
            pn = Yi[Sn | x] ^ nr;
            vn = Yi[Sn | $i] ^ sr;
            dn = Yi[Sn | tn] ^ ar;
            gn = Yi[Sn | en] ^ hr;
            bn = Yi[Sn | rn] ^ or;
            mn = Yi[Sn | nn] ^ fr;
            t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
            e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
            r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
            m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
            S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
            A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
            E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
            _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
            k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
            U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
            x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
            $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
            tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
            en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
            rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
            nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
            sn = Yi[Sn | t] ^ Me;
            an = Yi[Sn | e] ^ Ce;
            hn = Yi[Sn | r] ^ Ze;
            on = Yi[Sn | m] ^ ze;
            fn = Yi[Sn | S] ^ Oe;
            un = Yi[Sn | A] ^ Te;
            ln = Yi[Sn | E] ^ Re;
            cn = Yi[Sn | _] ^ Be;
            wn = Yi[Sn | k] ^ Pe;
            yn = Yi[Sn | U] ^ Ke;
            pn = Yi[Sn | x] ^ Ne;
            vn = Yi[Sn | $i] ^ je;
            dn = Yi[Sn | tn] ^ De;
            gn = Yi[Sn | en] ^ Ge;
            bn = Yi[Sn | rn] ^ Fe;
            mn = Yi[Sn | nn] ^ We;
            t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
            e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
            r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
            m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
            S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
            A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
            E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
            _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
            k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
            U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
            x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
            $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
            tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
            en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
            rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
            nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
            sn = Yi[Sn | t] ^ ve;
            an = Yi[Sn | e] ^ de;
            hn = Yi[Sn | r] ^ ge;
            on = Yi[Sn | m] ^ be;
            fn = Yi[Sn | S] ^ me;
            un = Yi[Sn | A] ^ Se;
            ln = Yi[Sn | E] ^ Ae;
            cn = Yi[Sn | _] ^ Ee;
            wn = Yi[Sn | k] ^ _e;
            yn = Yi[Sn | U] ^ ke;
            pn = Yi[Sn | x] ^ Ue;
            vn = Yi[Sn | $i] ^ xe;
            dn = Yi[Sn | tn] ^ Le;
            gn = Yi[Sn | en] ^ He;
            bn = Yi[Sn | rn] ^ Ie;
            mn = Yi[Sn | nn] ^ qe;
            t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
            e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
            r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
            m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
            S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
            A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
            E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
            _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
            k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
            U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
            x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
            $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
            tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
            en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
            rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
            nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
            sn = Yi[Sn | t] ^ te;
            an = Yi[Sn | e] ^ ee;
            hn = Yi[Sn | r] ^ re;
            on = Yi[Sn | m] ^ ie;
            fn = Yi[Sn | S] ^ ne;
            un = Yi[Sn | A] ^ se;
            ln = Yi[Sn | E] ^ ae;
            cn = Yi[Sn | _] ^ he;
            wn = Yi[Sn | k] ^ oe;
            yn = Yi[Sn | U] ^ fe;
            pn = Yi[Sn | x] ^ ue;
            vn = Yi[Sn | $i] ^ le;
            dn = Yi[Sn | tn] ^ ce;
            gn = Yi[Sn | en] ^ we;
            bn = Yi[Sn | rn] ^ ye;
            mn = Yi[Sn | nn] ^ pe;
            t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
            e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
            r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
            m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
            S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
            A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
            E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
            _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
            k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
            U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
            x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
            $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
            tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
            en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
            rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
            nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
            sn = Yi[Sn | t] ^ Rt;
            an = Yi[Sn | e] ^ Bt;
            hn = Yi[Sn | r] ^ Pt;
            on = Yi[Sn | m] ^ Kt;
            fn = Yi[Sn | S] ^ Nt;
            un = Yi[Sn | A] ^ jt;
            ln = Yi[Sn | E] ^ Dt;
            cn = Yi[Sn | _] ^ Gt;
            wn = Yi[Sn | k] ^ Ft;
            yn = Yi[Sn | U] ^ Wt;
            pn = Yi[Sn | x] ^ Vt;
            vn = Yi[Sn | $i] ^ Jt;
            dn = Yi[Sn | tn] ^ Qt;
            gn = Yi[Sn | en] ^ Xt;
            bn = Yi[Sn | rn] ^ Yt;
            mn = Yi[Sn | nn] ^ $t;
            t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
            e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
            r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
            m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
            S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
            A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
            E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
            _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
            k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
            U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
            x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
            $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
            tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
            en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
            rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
            nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
            sn = Yi[Sn | t] ^ At;
            an = Yi[Sn | e] ^ Et;
            hn = Yi[Sn | r] ^ _t;
            on = Yi[Sn | m] ^ kt;
            fn = Yi[Sn | S] ^ Ut;
            un = Yi[Sn | A] ^ xt;
            ln = Yi[Sn | E] ^ Lt;
            cn = Yi[Sn | _] ^ Ht;
            wn = Yi[Sn | k] ^ It;
            yn = Yi[Sn | U] ^ qt;
            pn = Yi[Sn | x] ^ Mt;
            vn = Yi[Sn | $i] ^ Ct;
            dn = Yi[Sn | tn] ^ Zt;
            gn = Yi[Sn | en] ^ zt;
            bn = Yi[Sn | rn] ^ Ot;
            mn = Yi[Sn | nn] ^ Tt;
            t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
            e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
            r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
            m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
            S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
            A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
            E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
            _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
            k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
            U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
            x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
            $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
            tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
            en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
            rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
            nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
            sn = Yi[Sn | t] ^ at;
            an = Yi[Sn | e] ^ ht;
            hn = Yi[Sn | r] ^ ot;
            on = Yi[Sn | m] ^ ft;
            fn = Yi[Sn | S] ^ ut;
            un = Yi[Sn | A] ^ lt;
            ln = Yi[Sn | E] ^ ct;
            cn = Yi[Sn | _] ^ wt;
            wn = Yi[Sn | k] ^ yt;
            yn = Yi[Sn | U] ^ pt;
            pn = Yi[Sn | x] ^ vt;
            vn = Yi[Sn | $i] ^ dt;
            dn = Yi[Sn | tn] ^ gt;
            gn = Yi[Sn | en] ^ bt;
            bn = Yi[Sn | rn] ^ mt;
            mn = Yi[Sn | nn] ^ St;
            t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
            e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
            r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
            m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
            S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
            A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
            E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
            _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
            k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
            U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
            x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
            $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
            tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
            en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
            rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
            nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
            sn = Yi[Sn | t] ^ D;
            an = Yi[Sn | e] ^ G;
            hn = Yi[Sn | r] ^ F;
            on = Yi[Sn | m] ^ W;
            fn = Yi[Sn | S] ^ V;
            un = Yi[Sn | A] ^ J;
            ln = Yi[Sn | E] ^ Q;
            cn = Yi[Sn | _] ^ X;
            wn = Yi[Sn | k] ^ Y;
            yn = Yi[Sn | U] ^ $;
            pn = Yi[Sn | x] ^ tt;
            vn = Yi[Sn | $i] ^ et;
            dn = Yi[Sn | tn] ^ rt;
            gn = Yi[Sn | en] ^ it;
            bn = Yi[Sn | rn] ^ nt;
            mn = Yi[Sn | nn] ^ st;
            t = Yi[kn | sn] ^ Yi[En | an] ^ Yi[_n | hn] ^ Yi[An | on];
            e = Yi[An | dn] ^ Yi[kn | gn] ^ Yi[En | bn] ^ Yi[_n | mn];
            r = Yi[_n | wn] ^ Yi[An | yn] ^ Yi[kn | pn] ^ Yi[En | vn];
            m = Yi[En | fn] ^ Yi[_n | un] ^ Yi[An | ln] ^ Yi[kn | cn];
            S = Yi[kn | fn] ^ Yi[En | un] ^ Yi[_n | ln] ^ Yi[An | cn];
            A = Yi[An | sn] ^ Yi[kn | an] ^ Yi[En | hn] ^ Yi[_n | on];
            E = Yi[_n | dn] ^ Yi[An | gn] ^ Yi[kn | bn] ^ Yi[En | mn];
            _ = Yi[En | wn] ^ Yi[_n | yn] ^ Yi[An | pn] ^ Yi[kn | vn];
            k = Yi[kn | wn] ^ Yi[En | yn] ^ Yi[_n | pn] ^ Yi[An | vn];
            U = Yi[An | fn] ^ Yi[kn | un] ^ Yi[En | ln] ^ Yi[_n | cn];
            x = Yi[_n | sn] ^ Yi[An | an] ^ Yi[kn | hn] ^ Yi[En | on];
            $i = Yi[En | dn] ^ Yi[_n | gn] ^ Yi[An | bn] ^ Yi[kn | mn];
            tn = Yi[kn | dn] ^ Yi[En | gn] ^ Yi[_n | bn] ^ Yi[An | mn];
            en = Yi[An | wn] ^ Yi[kn | yn] ^ Yi[En | pn] ^ Yi[_n | vn];
            rn = Yi[_n | fn] ^ Yi[An | un] ^ Yi[kn | ln] ^ Yi[En | cn];
            nn = Yi[En | sn] ^ Yi[_n | an] ^ Yi[An | hn] ^ Yi[kn | on];
            i = Yi[Sn | t] ^ L;
            n = Yi[Sn | e] ^ H;
            s = Yi[Sn | r] ^ I;
            a = Yi[Sn | m] ^ q;
            h = Yi[Sn | S] ^ M;
            o = Yi[Sn | A] ^ C;
            f = Yi[Sn | E] ^ Z;
            u = Yi[Sn | _] ^ z;
            l = Yi[Sn | k] ^ O;
            c = Yi[Sn | U] ^ T;
            w = Yi[Sn | x] ^ R;
            y = Yi[Sn | $i] ^ B;
            p = Yi[Sn | tn] ^ P;
            v = Yi[Sn | en] ^ K;
            d = Yi[Sn | rn] ^ N;
            g = Yi[Sn | nn] ^ j;
        }
        function nn(t, e, r, b, m, S, A, E, _, k, U, x, L, H, I, q) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            b = b | 0;
            m = m | 0;
            S = S | 0;
            A = A | 0;
            E = E | 0;
            _ = _ | 0;
            k = k | 0;
            U = U | 0;
            x = x | 0;
            L = L | 0;
            H = H | 0;
            I = I | 0;
            q = q | 0;
            i = t;
            n = e;
            s = r;
            a = b;
            h = m;
            o = S;
            f = A;
            u = E;
            l = _;
            c = k;
            w = U;
            y = x;
            p = L;
            v = H;
            d = I;
            g = q;
        }
        function sn(t) {
            t = t | 0;
            Yi[t] = i;
            Yi[t | 1] = n;
            Yi[t | 2] = s;
            Yi[t | 3] = a;
            Yi[t | 4] = h;
            Yi[t | 5] = o;
            Yi[t | 6] = f;
            Yi[t | 7] = u;
            Yi[t | 8] = l;
            Yi[t | 9] = c;
            Yi[t | 10] = w;
            Yi[t | 11] = y;
            Yi[t | 12] = p;
            Yi[t | 13] = v;
            Yi[t | 14] = d;
            Yi[t | 15] = g;
        }
        function an(t, e, r, i, n, s, a, h, o, f, u, l, c, w, y, p) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            n = n | 0;
            s = s | 0;
            a = a | 0;
            h = h | 0;
            o = o | 0;
            f = f | 0;
            u = u | 0;
            l = l | 0;
            c = c | 0;
            w = w | 0;
            y = y | 0;
            p = p | 0;
            L = t;
            H = e;
            I = r;
            q = i;
            M = n;
            C = s;
            Z = a;
            z = h;
            O = o;
            T = f;
            R = u;
            B = l;
            P = c;
            K = w;
            N = y;
            j = p;
            b = 16;
            $i();
        }
        function hn(t, e, r, i, n, s, a, h, o, f, u, l, c, w, y, p, v, d, g, m, S, A, E, _, k, U, x, at, ht, ot, ft, ut) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            n = n | 0;
            s = s | 0;
            a = a | 0;
            h = h | 0;
            o = o | 0;
            f = f | 0;
            u = u | 0;
            l = l | 0;
            c = c | 0;
            w = w | 0;
            y = y | 0;
            p = p | 0;
            v = v | 0;
            d = d | 0;
            g = g | 0;
            m = m | 0;
            S = S | 0;
            A = A | 0;
            E = E | 0;
            _ = _ | 0;
            k = k | 0;
            U = U | 0;
            x = x | 0;
            at = at | 0;
            ht = ht | 0;
            ot = ot | 0;
            ft = ft | 0;
            ut = ut | 0;
            L = t;
            H = e;
            I = r;
            q = i;
            M = n;
            C = s;
            Z = a;
            z = h;
            O = o;
            T = f;
            R = u;
            B = l;
            P = c;
            K = w;
            N = y;
            j = p;
            D = v;
            G = d;
            F = g;
            W = m;
            V = S;
            J = A;
            Q = E;
            X = _;
            Y = k;
            $ = U;
            tt = x;
            et = at;
            rt = ht;
            it = ot;
            nt = ft;
            st = ut;
            b = 32;
            tn();
        }
        function on(t, e) {
            t = t | 0;
            e = e | 0;
            var r = 0;
            if (t & 15) return -1;
            while ((e | 0) >= 16) {
                en(Yi[t] | 0, Yi[t | 1] | 0, Yi[t | 2] | 0, Yi[t | 3] | 0, Yi[t | 4] | 0, Yi[t | 5] | 0, Yi[t | 6] | 0, Yi[t | 7] | 0, Yi[t | 8] | 0, Yi[t | 9] | 0, Yi[t | 10] | 0, Yi[t | 11] | 0, Yi[t | 12] | 0, Yi[t | 13] | 0, Yi[t | 14] | 0, Yi[t | 15] | 0);
                Yi[t] = i;
                Yi[t | 1] = n;
                Yi[t | 2] = s;
                Yi[t | 3] = a;
                Yi[t | 4] = h;
                Yi[t | 5] = o;
                Yi[t | 6] = f;
                Yi[t | 7] = u;
                Yi[t | 8] = l;
                Yi[t | 9] = c;
                Yi[t | 10] = w;
                Yi[t | 11] = y;
                Yi[t | 12] = p;
                Yi[t | 13] = v;
                Yi[t | 14] = d;
                Yi[t | 15] = g;
                t = t + 16 | 0;
                e = e - 16 | 0;
                r = r + 16 | 0;
            }
            return r | 0;
        }
        function fn(t, e) {
            t = t | 0;
            e = e | 0;
            var r = 0;
            if (t & 15) return -1;
            while ((e | 0) >= 16) {
                rn(Yi[t] | 0, Yi[t | 1] | 0, Yi[t | 2] | 0, Yi[t | 3] | 0, Yi[t | 4] | 0, Yi[t | 5] | 0, Yi[t | 6] | 0, Yi[t | 7] | 0, Yi[t | 8] | 0, Yi[t | 9] | 0, Yi[t | 10] | 0, Yi[t | 11] | 0, Yi[t | 12] | 0, Yi[t | 13] | 0, Yi[t | 14] | 0, Yi[t | 15] | 0);
                Yi[t] = i;
                Yi[t | 1] = n;
                Yi[t | 2] = s;
                Yi[t | 3] = a;
                Yi[t | 4] = h;
                Yi[t | 5] = o;
                Yi[t | 6] = f;
                Yi[t | 7] = u;
                Yi[t | 8] = l;
                Yi[t | 9] = c;
                Yi[t | 10] = w;
                Yi[t | 11] = y;
                Yi[t | 12] = p;
                Yi[t | 13] = v;
                Yi[t | 14] = d;
                Yi[t | 15] = g;
                t = t + 16 | 0;
                e = e - 16 | 0;
                r = r + 16 | 0;
            }
            return r | 0;
        }
        function un(t, e) {
            t = t | 0;
            e = e | 0;
            var r = 0;
            if (t & 15) return -1;
            while ((e | 0) >= 16) {
                en(i ^ Yi[t], n ^ Yi[t | 1], s ^ Yi[t | 2], a ^ Yi[t | 3], h ^ Yi[t | 4], o ^ Yi[t | 5], f ^ Yi[t | 6], u ^ Yi[t | 7], l ^ Yi[t | 8], c ^ Yi[t | 9], w ^ Yi[t | 10], y ^ Yi[t | 11], p ^ Yi[t | 12], v ^ Yi[t | 13], d ^ Yi[t | 14], g ^ Yi[t | 15]);
                Yi[t] = i;
                Yi[t | 1] = n;
                Yi[t | 2] = s;
                Yi[t | 3] = a;
                Yi[t | 4] = h;
                Yi[t | 5] = o;
                Yi[t | 6] = f;
                Yi[t | 7] = u;
                Yi[t | 8] = l;
                Yi[t | 9] = c;
                Yi[t | 10] = w;
                Yi[t | 11] = y;
                Yi[t | 12] = p;
                Yi[t | 13] = v;
                Yi[t | 14] = d;
                Yi[t | 15] = g;
                t = t + 16 | 0;
                e = e - 16 | 0;
                r = r + 16 | 0;
            }
            return r | 0;
        }
        function ln(t, e) {
            t = t | 0;
            e = e | 0;
            var r = 0, b = 0, m = 0, S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0, L = 0, H = 0, I = 0, q = 0, M = 0, C = 0, Z = 0;
            if (t & 15) return -1;
            r = i;
            b = n;
            m = s;
            S = a;
            A = h;
            E = o;
            _ = f;
            k = u;
            U = l;
            x = c;
            L = w;
            H = y;
            I = p;
            q = v;
            M = d;
            C = g;
            while ((e | 0) >= 16) {
                rn(Yi[t] | 0, Yi[t | 1] | 0, Yi[t | 2] | 0, Yi[t | 3] | 0, Yi[t | 4] | 0, Yi[t | 5] | 0, Yi[t | 6] | 0, Yi[t | 7] | 0, Yi[t | 8] | 0, Yi[t | 9] | 0, Yi[t | 10] | 0, Yi[t | 11] | 0, Yi[t | 12] | 0, Yi[t | 13] | 0, Yi[t | 14] | 0, Yi[t | 15] | 0);
                i = i ^ r;
                r = Yi[t] | 0;
                n = n ^ b;
                b = Yi[t | 1] | 0;
                s = s ^ m;
                m = Yi[t | 2] | 0;
                a = a ^ S;
                S = Yi[t | 3] | 0;
                h = h ^ A;
                A = Yi[t | 4] | 0;
                o = o ^ E;
                E = Yi[t | 5] | 0;
                f = f ^ _;
                _ = Yi[t | 6] | 0;
                u = u ^ k;
                k = Yi[t | 7] | 0;
                l = l ^ U;
                U = Yi[t | 8] | 0;
                c = c ^ x;
                x = Yi[t | 9] | 0;
                w = w ^ L;
                L = Yi[t | 10] | 0;
                y = y ^ H;
                H = Yi[t | 11] | 0;
                p = p ^ I;
                I = Yi[t | 12] | 0;
                v = v ^ q;
                q = Yi[t | 13] | 0;
                d = d ^ M;
                M = Yi[t | 14] | 0;
                g = g ^ C;
                C = Yi[t | 15] | 0;
                Yi[t] = i;
                Yi[t | 1] = n;
                Yi[t | 2] = s;
                Yi[t | 3] = a;
                Yi[t | 4] = h;
                Yi[t | 5] = o;
                Yi[t | 6] = f;
                Yi[t | 7] = u;
                Yi[t | 8] = l;
                Yi[t | 9] = c;
                Yi[t | 10] = w;
                Yi[t | 11] = y;
                Yi[t | 12] = p;
                Yi[t | 13] = v;
                Yi[t | 14] = d;
                Yi[t | 15] = g;
                t = t + 16 | 0;
                e = e - 16 | 0;
                Z = Z + 16 | 0;
            }
            i = r;
            n = b;
            s = m;
            a = S;
            h = A;
            o = E;
            f = _;
            u = k;
            l = U;
            c = x;
            w = L;
            y = H;
            p = I;
            v = q;
            d = M;
            g = C;
            return Z | 0;
        }
        function cn(t, e, r) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            if (t & 15) return -1;
            if (~r) if (r & 31) return -1;
            while ((e | 0) >= 16) {
                en(i ^ Yi[t], n ^ Yi[t | 1], s ^ Yi[t | 2], a ^ Yi[t | 3], h ^ Yi[t | 4], o ^ Yi[t | 5], f ^ Yi[t | 6], u ^ Yi[t | 7], l ^ Yi[t | 8], c ^ Yi[t | 9], w ^ Yi[t | 10], y ^ Yi[t | 11], p ^ Yi[t | 12], v ^ Yi[t | 13], d ^ Yi[t | 14], g ^ Yi[t | 15]);
                t = t + 16 | 0;
                e = e - 16 | 0;
            }
            if ((e | 0) > 0) {
                i = i ^ Yi[t];
                if ((e | 0) > 1) n = n ^ Yi[t | 1];
                if ((e | 0) > 2) s = s ^ Yi[t | 2];
                if ((e | 0) > 3) a = a ^ Yi[t | 3];
                if ((e | 0) > 4) h = h ^ Yi[t | 4];
                if ((e | 0) > 5) o = o ^ Yi[t | 5];
                if ((e | 0) > 6) f = f ^ Yi[t | 6];
                if ((e | 0) > 7) u = u ^ Yi[t | 7];
                if ((e | 0) > 8) l = l ^ Yi[t | 8];
                if ((e | 0) > 9) c = c ^ Yi[t | 9];
                if ((e | 0) > 10) w = w ^ Yi[t | 10];
                if ((e | 0) > 11) y = y ^ Yi[t | 11];
                if ((e | 0) > 12) p = p ^ Yi[t | 12];
                if ((e | 0) > 13) v = v ^ Yi[t | 13];
                if ((e | 0) > 14) d = d ^ Yi[t | 14];
                en(i, n, s, a, h, o, f, u, l, c, w, y, p, v, d, g);
                t = t + e | 0;
                e = 0;
            }
            if (~r) {
                Yi[r | 0] = i;
                Yi[r | 1] = n;
                Yi[r | 2] = s;
                Yi[r | 3] = a;
                Yi[r | 4] = h;
                Yi[r | 5] = o;
                Yi[r | 6] = f;
                Yi[r | 7] = u;
                Yi[r | 8] = l;
                Yi[r | 9] = c;
                Yi[r | 10] = w;
                Yi[r | 11] = y;
                Yi[r | 12] = p;
                Yi[r | 13] = v;
                Yi[r | 14] = d;
                Yi[r | 15] = g;
            }
            return 0;
        }
        function wn(t, e, r, b, m, S, A, E, _, k, U, x, L, H, I) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            b = b | 0;
            m = m | 0;
            S = S | 0;
            A = A | 0;
            E = E | 0;
            _ = _ | 0;
            k = k | 0;
            U = U | 0;
            x = x | 0;
            L = L | 0;
            H = H | 0;
            I = I | 0;
            var q = 0;
            while ((e | 0) >= 16) {
                en(r, b, m, S, A, E, _, k, U, x, L, H, I >>> 24, I >>> 16 & 255, I >>> 8 & 255, I & 255);
                Yi[t | 0] = Yi[t | 0] ^ i;
                Yi[t | 1] = Yi[t | 1] ^ n;
                Yi[t | 2] = Yi[t | 2] ^ s;
                Yi[t | 3] = Yi[t | 3] ^ a;
                Yi[t | 4] = Yi[t | 4] ^ h;
                Yi[t | 5] = Yi[t | 5] ^ o;
                Yi[t | 6] = Yi[t | 6] ^ f;
                Yi[t | 7] = Yi[t | 7] ^ u;
                Yi[t | 8] = Yi[t | 8] ^ l;
                Yi[t | 9] = Yi[t | 9] ^ c;
                Yi[t | 10] = Yi[t | 10] ^ w;
                Yi[t | 11] = Yi[t | 11] ^ y;
                Yi[t | 12] = Yi[t | 12] ^ p;
                Yi[t | 13] = Yi[t | 13] ^ v;
                Yi[t | 14] = Yi[t | 14] ^ d;
                Yi[t | 15] = Yi[t | 15] ^ g;
                t = t + 16 | 0;
                e = e - 16 | 0;
                q = q + 16 | 0;
                I = I + 1 | 0;
            }
            return q | 0;
        }
        function yn(t, e, r, b, m, S, A, E, _, k, U, x, L, H, I, q, M, C) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            b = b | 0;
            m = m | 0;
            S = S | 0;
            A = A | 0;
            E = E | 0;
            _ = _ | 0;
            k = k | 0;
            U = U | 0;
            x = x | 0;
            L = L | 0;
            H = H | 0;
            I = I | 0;
            q = q | 0;
            M = M | 0;
            C = C | 0;
            var Z = 0, z = 0, O = 0, T = 0, R = 0, B = 0, P = 0, K = 0, N = 0, j = 0, D = 0, G = 0, F = 0, W = 0, V = 0, J = 0, Q = 0, X = 0, Y = 0, $ = 0, tt = 0, et = 0, rt = 0, it = 0, nt = 0, st = 0, at = 0, ht = 0, ot = 0, ft = 0, ut = 0, lt = 0, ct = 0;
            if (t & 15) return -1;
            Z = i, z = n, O = s, T = a, R = h, B = o, P = f, K = u, N = l, j = c, D = w, G = y, 
            F = p, W = v, V = d, J = g;
            while ((e | 0) >= 16) {
                Q = Yi[t] | 0;
                X = Yi[t | 1] | 0;
                Y = Yi[t | 2] | 0;
                $ = Yi[t | 3] | 0;
                tt = Yi[t | 4] | 0;
                et = Yi[t | 5] | 0;
                rt = Yi[t | 6] | 0;
                it = Yi[t | 7] | 0;
                nt = Yi[t | 8] | 0;
                st = Yi[t | 9] | 0;
                at = Yi[t | 10] | 0;
                ht = Yi[t | 11] | 0;
                ot = Yi[t | 12] | 0;
                ft = Yi[t | 13] | 0;
                ut = Yi[t | 14] | 0;
                lt = Yi[t | 15] | 0;
                en(r, b, m, S, A, E, _, k, U ^ M >>> 24, x ^ M >>> 16 & 255, L ^ M >>> 8 & 255, H ^ M & 255, I ^ C >>> 24, q ^ C >>> 16 & 255, C >>> 8 & 255, C & 255);
                Yi[t] = Q ^ i;
                Yi[t | 1] = X ^ n;
                Yi[t | 2] = Y ^ s;
                Yi[t | 3] = $ ^ a;
                Yi[t | 4] = tt ^ h;
                Yi[t | 5] = et ^ o;
                Yi[t | 6] = rt ^ f;
                Yi[t | 7] = it ^ u;
                Yi[t | 8] = nt ^ l;
                Yi[t | 9] = st ^ c;
                Yi[t | 10] = at ^ w;
                Yi[t | 11] = ht ^ y;
                Yi[t | 12] = ot ^ p;
                Yi[t | 13] = ft ^ v;
                Yi[t | 14] = ut ^ d;
                Yi[t | 15] = lt ^ g;
                en(Q ^ Z, X ^ z, Y ^ O, $ ^ T, tt ^ R, et ^ B, rt ^ P, it ^ K, nt ^ N, st ^ j, at ^ D, ht ^ G, ot ^ F, ft ^ W, ut ^ V, lt ^ J);
                Z = i, z = n, O = s, T = a, R = h, B = o, P = f, K = u, N = l, j = c, D = w, G = y, 
                F = p, W = v, V = d, J = g;
                ct = ct + 16 | 0;
                t = t + 16 | 0;
                e = e - 16 | 0;
                C = C + 1 | 0;
                if ((C | 0) == 0) M = M + 1 | 0;
            }
            if ((e | 0) > 0) {
                Q = Yi[t] | 0;
                X = (e | 0) > 1 ? Yi[t | 1] | 0 : 0;
                Y = (e | 0) > 2 ? Yi[t | 2] | 0 : 0;
                $ = (e | 0) > 3 ? Yi[t | 3] | 0 : 0;
                tt = (e | 0) > 4 ? Yi[t | 4] | 0 : 0;
                et = (e | 0) > 5 ? Yi[t | 5] | 0 : 0;
                rt = (e | 0) > 6 ? Yi[t | 6] | 0 : 0;
                it = (e | 0) > 7 ? Yi[t | 7] | 0 : 0;
                nt = (e | 0) > 8 ? Yi[t | 8] | 0 : 0;
                st = (e | 0) > 9 ? Yi[t | 9] | 0 : 0;
                at = (e | 0) > 10 ? Yi[t | 10] | 0 : 0;
                ht = (e | 0) > 11 ? Yi[t | 11] | 0 : 0;
                ot = (e | 0) > 12 ? Yi[t | 12] | 0 : 0;
                ft = (e | 0) > 13 ? Yi[t | 13] | 0 : 0;
                ut = (e | 0) > 14 ? Yi[t | 14] | 0 : 0;
                en(r, b, m, S, A, E, _, k, U ^ M >>> 24, x ^ M >>> 16 & 255, L ^ M >>> 8 & 255, H ^ M & 255, I ^ C >>> 24, q ^ C >>> 16 & 255, C >>> 8 & 255, C & 255);
                Yi[t] = Q ^ i;
                if ((e | 0) > 1) Yi[t | 1] = X ^ n;
                if ((e | 0) > 2) Yi[t | 2] = Y ^ s;
                if ((e | 0) > 3) Yi[t | 3] = $ ^ a;
                if ((e | 0) > 4) Yi[t | 4] = tt ^ h;
                if ((e | 0) > 5) Yi[t | 5] = et ^ o;
                if ((e | 0) > 6) Yi[t | 6] = rt ^ f;
                if ((e | 0) > 7) Yi[t | 7] = it ^ u;
                if ((e | 0) > 8) Yi[t | 8] = nt ^ l;
                if ((e | 0) > 9) Yi[t | 9] = st ^ c;
                if ((e | 0) > 10) Yi[t | 10] = at ^ w;
                if ((e | 0) > 11) Yi[t | 11] = ht ^ y;
                if ((e | 0) > 12) Yi[t | 12] = ot ^ p;
                if ((e | 0) > 13) Yi[t | 13] = ft ^ v;
                if ((e | 0) > 14) Yi[t | 14] = ut ^ d;
                en(Q ^ Z, X ^ z, Y ^ O, $ ^ T, tt ^ R, et ^ B, rt ^ P, it ^ K, nt ^ N, st ^ j, at ^ D, ht ^ G, ot ^ F, ft ^ W, ut ^ V, J);
                Z = i, z = n, O = s, T = a, R = h, B = o, P = f, K = u, N = l, j = c, D = w, G = y, 
                F = p, W = v, V = d, J = g;
                ct = ct + e | 0;
                t = t + e | 0;
                e = 0;
                C = C + 1 | 0;
                if ((C | 0) == 0) M = M + 1 | 0;
            }
            return ct | 0;
        }
        function pn(t, e, r, b, m, S, A, E, _, k, U, x, L, H, I, q, M, C) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            b = b | 0;
            m = m | 0;
            S = S | 0;
            A = A | 0;
            E = E | 0;
            _ = _ | 0;
            k = k | 0;
            U = U | 0;
            x = x | 0;
            L = L | 0;
            H = H | 0;
            I = I | 0;
            q = q | 0;
            M = M | 0;
            C = C | 0;
            var Z = 0, z = 0, O = 0, T = 0, R = 0, B = 0, P = 0, K = 0, N = 0, j = 0, D = 0, G = 0, F = 0, W = 0, V = 0, J = 0, Q = 0, X = 0, Y = 0, $ = 0, tt = 0, et = 0, rt = 0, it = 0, nt = 0, st = 0, at = 0, ht = 0, ot = 0, ft = 0, ut = 0, lt = 0, ct = 0;
            if (t & 15) return -1;
            Z = i, z = n, O = s, T = a, R = h, B = o, P = f, K = u, N = l, j = c, D = w, G = y, 
            F = p, W = v, V = d, J = g;
            while ((e | 0) >= 16) {
                en(r, b, m, S, A, E, _, k, U ^ M >>> 24, x ^ M >>> 16 & 255, L ^ M >>> 8 & 255, H ^ M & 255, I ^ C >>> 24, q ^ C >>> 16 & 255, C >>> 8 & 255, C & 255);
                Yi[t] = Q = Yi[t] ^ i;
                Yi[t | 1] = X = Yi[t | 1] ^ n;
                Yi[t | 2] = Y = Yi[t | 2] ^ s;
                Yi[t | 3] = $ = Yi[t | 3] ^ a;
                Yi[t | 4] = tt = Yi[t | 4] ^ h;
                Yi[t | 5] = et = Yi[t | 5] ^ o;
                Yi[t | 6] = rt = Yi[t | 6] ^ f;
                Yi[t | 7] = it = Yi[t | 7] ^ u;
                Yi[t | 8] = nt = Yi[t | 8] ^ l;
                Yi[t | 9] = st = Yi[t | 9] ^ c;
                Yi[t | 10] = at = Yi[t | 10] ^ w;
                Yi[t | 11] = ht = Yi[t | 11] ^ y;
                Yi[t | 12] = ot = Yi[t | 12] ^ p;
                Yi[t | 13] = ft = Yi[t | 13] ^ v;
                Yi[t | 14] = ut = Yi[t | 14] ^ d;
                Yi[t | 15] = lt = Yi[t | 15] ^ g;
                en(Q ^ Z, X ^ z, Y ^ O, $ ^ T, tt ^ R, et ^ B, rt ^ P, it ^ K, nt ^ N, st ^ j, at ^ D, ht ^ G, ot ^ F, ft ^ W, ut ^ V, lt ^ J);
                Z = i, z = n, O = s, T = a, R = h, B = o, P = f, K = u, N = l, j = c, D = w, G = y, 
                F = p, W = v, V = d, J = g;
                ct = ct + 16 | 0;
                t = t + 16 | 0;
                e = e - 16 | 0;
                C = C + 1 | 0;
                if ((C | 0) == 0) M = M + 1 | 0;
            }
            if ((e | 0) > 0) {
                en(r, b, m, S, A, E, _, k, U ^ M >>> 24, x ^ M >>> 16 & 255, L ^ M >>> 8 & 255, H ^ M & 255, I ^ C >>> 24, q ^ C >>> 16 & 255, C >>> 8 & 255, C & 255);
                Q = Yi[t] ^ i;
                X = (e | 0) > 1 ? Yi[t | 1] ^ n : 0;
                Y = (e | 0) > 2 ? Yi[t | 2] ^ s : 0;
                $ = (e | 0) > 3 ? Yi[t | 3] ^ a : 0;
                tt = (e | 0) > 4 ? Yi[t | 4] ^ h : 0;
                et = (e | 0) > 5 ? Yi[t | 5] ^ o : 0;
                rt = (e | 0) > 6 ? Yi[t | 6] ^ f : 0;
                it = (e | 0) > 7 ? Yi[t | 7] ^ u : 0;
                nt = (e | 0) > 8 ? Yi[t | 8] ^ l : 0;
                st = (e | 0) > 9 ? Yi[t | 9] ^ c : 0;
                at = (e | 0) > 10 ? Yi[t | 10] ^ w : 0;
                ht = (e | 0) > 11 ? Yi[t | 11] ^ y : 0;
                ot = (e | 0) > 12 ? Yi[t | 12] ^ p : 0;
                ft = (e | 0) > 13 ? Yi[t | 13] ^ v : 0;
                ut = (e | 0) > 14 ? Yi[t | 14] ^ d : 0;
                lt = (e | 0) > 15 ? Yi[t | 15] ^ g : 0;
                Yi[t] = Q;
                if ((e | 0) > 1) Yi[t | 1] = X;
                if ((e | 0) > 2) Yi[t | 2] = Y;
                if ((e | 0) > 3) Yi[t | 3] = $;
                if ((e | 0) > 4) Yi[t | 4] = tt;
                if ((e | 0) > 5) Yi[t | 5] = et;
                if ((e | 0) > 6) Yi[t | 6] = rt;
                if ((e | 0) > 7) Yi[t | 7] = it;
                if ((e | 0) > 8) Yi[t | 8] = nt;
                if ((e | 0) > 9) Yi[t | 9] = st;
                if ((e | 0) > 10) Yi[t | 10] = at;
                if ((e | 0) > 11) Yi[t | 11] = ht;
                if ((e | 0) > 12) Yi[t | 12] = ot;
                if ((e | 0) > 13) Yi[t | 13] = ft;
                if ((e | 0) > 14) Yi[t | 14] = ut;
                en(Q ^ Z, X ^ z, Y ^ O, $ ^ T, tt ^ R, et ^ B, rt ^ P, it ^ K, nt ^ N, st ^ j, at ^ D, ht ^ G, ot ^ F, ft ^ W, ut ^ V, lt ^ J);
                Z = i, z = n, O = s, T = a, R = h, B = o, P = f, K = u, N = l, j = c, D = w, G = y, 
                F = p, W = v, V = d, J = g;
                ct = ct + e | 0;
                t = t + e | 0;
                e = 0;
                C = C + 1 | 0;
                if ((C | 0) == 0) M = M + 1 | 0;
            }
            return ct | 0;
        }
        function vn(t, e) {
            t = t | 0;
            e = e | 0;
            var r = 0;
            if (t & 15) return -1;
            while ((e | 0) >= 16) {
                en(i, n, s, a, h, o, f, u, l, c, w, y, p, v, d, g);
                i = i ^ Yi[t];
                n = n ^ Yi[t | 1];
                s = s ^ Yi[t | 2];
                a = a ^ Yi[t | 3];
                h = h ^ Yi[t | 4];
                o = o ^ Yi[t | 5];
                f = f ^ Yi[t | 6];
                u = u ^ Yi[t | 7];
                l = l ^ Yi[t | 8];
                c = c ^ Yi[t | 9];
                w = w ^ Yi[t | 10];
                y = y ^ Yi[t | 11];
                p = p ^ Yi[t | 12];
                v = v ^ Yi[t | 13];
                d = d ^ Yi[t | 14];
                g = g ^ Yi[t | 15];
                Yi[t] = i;
                Yi[t | 1] = n;
                Yi[t | 2] = s;
                Yi[t | 3] = a;
                Yi[t | 4] = h;
                Yi[t | 5] = o;
                Yi[t | 6] = f;
                Yi[t | 7] = u;
                Yi[t | 8] = l;
                Yi[t | 9] = c;
                Yi[t | 10] = w;
                Yi[t | 11] = y;
                Yi[t | 12] = p;
                Yi[t | 13] = v;
                Yi[t | 14] = d;
                Yi[t | 15] = g;
                t = t + 16 | 0;
                e = e - 16 | 0;
                r = r + 16 | 0;
            }
            if ((e | 0) > 0) {
                en(i, n, s, a, h, o, f, u, l, c, w, y, p, v, d, g);
                Yi[t] = Yi[t] ^ i;
                if ((e | 0) > 1) Yi[t | 1] = Yi[t | 1] ^ n;
                if ((e | 0) > 2) Yi[t | 2] = Yi[t | 2] ^ s;
                if ((e | 0) > 3) Yi[t | 3] = Yi[t | 3] ^ a;
                if ((e | 0) > 4) Yi[t | 4] = Yi[t | 4] ^ h;
                if ((e | 0) > 5) Yi[t | 5] = Yi[t | 5] ^ o;
                if ((e | 0) > 6) Yi[t | 6] = Yi[t | 6] ^ f;
                if ((e | 0) > 7) Yi[t | 7] = Yi[t | 7] ^ u;
                if ((e | 0) > 8) Yi[t | 8] = Yi[t | 8] ^ l;
                if ((e | 0) > 9) Yi[t | 9] = Yi[t | 9] ^ c;
                if ((e | 0) > 10) Yi[t | 10] = Yi[t | 10] ^ w;
                if ((e | 0) > 11) Yi[t | 11] = Yi[t | 11] ^ y;
                if ((e | 0) > 12) Yi[t | 12] = Yi[t | 12] ^ p;
                if ((e | 0) > 13) Yi[t | 13] = Yi[t | 13] ^ v;
                if ((e | 0) > 14) Yi[t | 14] = Yi[t | 14] ^ d;
                r = r + e | 0;
                t = t + e | 0;
                e = 0;
            }
            return r | 0;
        }
        function dn(t, e) {
            t = t | 0;
            e = e | 0;
            var r = 0, b = 0, m = 0, S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0, L = 0, H = 0, I = 0, q = 0, M = 0, C = 0, Z = 0;
            if (t & 15) return -1;
            while ((e | 0) >= 16) {
                en(i, n, s, a, h, o, f, u, l, c, w, y, p, v, d, g);
                r = Yi[t] | 0;
                b = Yi[t | 1] | 0;
                m = Yi[t | 2] | 0;
                S = Yi[t | 3] | 0;
                A = Yi[t | 4] | 0;
                E = Yi[t | 5] | 0;
                _ = Yi[t | 6] | 0;
                k = Yi[t | 7] | 0;
                U = Yi[t | 8] | 0;
                x = Yi[t | 9] | 0;
                L = Yi[t | 10] | 0;
                H = Yi[t | 11] | 0;
                I = Yi[t | 12] | 0;
                q = Yi[t | 13] | 0;
                M = Yi[t | 14] | 0;
                C = Yi[t | 15] | 0;
                Yi[t] = i ^ r;
                Yi[t | 1] = n ^ b;
                Yi[t | 2] = s ^ m;
                Yi[t | 3] = a ^ S;
                Yi[t | 4] = h ^ A;
                Yi[t | 5] = o ^ E;
                Yi[t | 6] = f ^ _;
                Yi[t | 7] = u ^ k;
                Yi[t | 8] = l ^ U;
                Yi[t | 9] = c ^ x;
                Yi[t | 10] = w ^ L;
                Yi[t | 11] = y ^ H;
                Yi[t | 12] = p ^ I;
                Yi[t | 13] = v ^ q;
                Yi[t | 14] = d ^ M;
                Yi[t | 15] = g ^ C;
                i = r;
                n = b;
                s = m;
                a = S;
                h = A;
                o = E;
                f = _;
                u = k;
                l = U;
                c = x;
                w = L;
                y = H;
                p = I;
                v = q;
                d = M;
                g = C;
                t = t + 16 | 0;
                e = e - 16 | 0;
                Z = Z + 16 | 0;
            }
            if ((e | 0) > 0) {
                en(i, n, s, a, h, o, f, u, l, c, w, y, p, v, d, g);
                Yi[t] = Yi[t] ^ i;
                if ((e | 0) > 1) Yi[t | 1] = Yi[t | 1] ^ n;
                if ((e | 0) > 2) Yi[t | 2] = Yi[t | 2] ^ s;
                if ((e | 0) > 3) Yi[t | 3] = Yi[t | 3] ^ a;
                if ((e | 0) > 4) Yi[t | 4] = Yi[t | 4] ^ h;
                if ((e | 0) > 5) Yi[t | 5] = Yi[t | 5] ^ o;
                if ((e | 0) > 6) Yi[t | 6] = Yi[t | 6] ^ f;
                if ((e | 0) > 7) Yi[t | 7] = Yi[t | 7] ^ u;
                if ((e | 0) > 8) Yi[t | 8] = Yi[t | 8] ^ l;
                if ((e | 0) > 9) Yi[t | 9] = Yi[t | 9] ^ c;
                if ((e | 0) > 10) Yi[t | 10] = Yi[t | 10] ^ w;
                if ((e | 0) > 11) Yi[t | 11] = Yi[t | 11] ^ y;
                if ((e | 0) > 12) Yi[t | 12] = Yi[t | 12] ^ p;
                if ((e | 0) > 13) Yi[t | 13] = Yi[t | 13] ^ v;
                if ((e | 0) > 14) Yi[t | 14] = Yi[t | 14] ^ d;
                Z = Z + e | 0;
                t = t + e | 0;
                e = 0;
            }
            return Z | 0;
        }
        function gn(t, e, r, i) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            var n = 0, s = 0, a = 0, h = 0, o = 0, f = 0, u = 0, l = 0, c = 0, w = 0;
            n = m | 0, s = S | 0, a = A | 0, h = E | 0;
            for (;(c | 0) < 128; c = c + 1 | 0) {
                if (n >>> 31) {
                    o = o ^ t, f = f ^ e, u = u ^ r, l = l ^ i;
                }
                n = n << 1 | s >>> 31, s = s << 1 | a >>> 31, a = a << 1 | h >>> 31, h = h << 1;
                w = i & 1;
                i = i >>> 1 | r << 31, r = r >>> 1 | e << 31, e = e >>> 1 | t << 31, t = t >>> 1;
                if (w) t = t ^ 3774873600;
            }
            _ = o, k = f, U = u, x = l;
        }
        function bn() {
            en(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), m = i << 24 | n << 16 | s << 8 | a, 
            S = h << 24 | o << 16 | f << 8 | u, A = l << 24 | c << 16 | w << 8 | y, E = p << 24 | v << 16 | d << 8 | g;
            _ = k = U = x = 0;
        }
        function mn(t, e) {
            t = t | 0;
            e = e | 0;
            var r = 0;
            if (t & 15) return -1;
            _ = i << 24 | n << 16 | s << 8 | a, k = h << 24 | o << 16 | f << 8 | u, U = l << 24 | c << 16 | w << 8 | y, 
            x = p << 24 | v << 16 | d << 8 | g;
            while ((e | 0) >= 16) {
                gn(_ ^ (Yi[t | 0] << 24 | Yi[t | 1] << 16 | Yi[t | 2] << 8 | Yi[t | 3]), k ^ (Yi[t | 4] << 24 | Yi[t | 5] << 16 | Yi[t | 6] << 8 | Yi[t | 7]), U ^ (Yi[t | 8] << 24 | Yi[t | 9] << 16 | Yi[t | 10] << 8 | Yi[t | 11]), x ^ (Yi[t | 12] << 24 | Yi[t | 13] << 16 | Yi[t | 14] << 8 | Yi[t | 15]));
                t = t + 16 | 0, e = e - 16 | 0, r = r + 16 | 0;
            }
            i = _ >>> 24, n = _ >>> 16 & 255, s = _ >>> 8 & 255, a = _ & 255, h = k >>> 24, 
            o = k >>> 16 & 255, f = k >>> 8 & 255, u = k & 255, l = U >>> 24, c = U >>> 16 & 255, 
            w = U >>> 8 & 255, y = U & 255, p = x >>> 24, v = x >>> 16 & 255, d = x >>> 8 & 255, 
            g = x & 255;
            return r | 0;
        }
        function Sn(t, e, r, b, m, S, A, E, L, H, I, q, M, C, Z) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            b = b | 0;
            m = m | 0;
            S = S | 0;
            A = A | 0;
            E = E | 0;
            L = L | 0;
            H = H | 0;
            I = I | 0;
            q = q | 0;
            M = M | 0;
            C = C | 0;
            Z = Z | 0;
            var z = 0, O = 0, T = 0, R = 0, B = 0, P = 0, K = 0, N = 0, j = 0, D = 0, G = 0, F = 0, W = 0, V = 0, J = 0, Q = 0, X = 0;
            if (t & 15) return -1;
            while ((e | 0) >= 16) {
                en(r, b, m, S, A, E, L, H, I, q, M, C, Z >>> 24, Z >>> 16 & 255, Z >>> 8 & 255, Z & 255);
                Yi[t | 0] = z = Yi[t | 0] ^ i, Yi[t | 1] = O = Yi[t | 1] ^ n, Yi[t | 2] = T = Yi[t | 2] ^ s, 
                Yi[t | 3] = R = Yi[t | 3] ^ a, Yi[t | 4] = B = Yi[t | 4] ^ h, Yi[t | 5] = P = Yi[t | 5] ^ o, 
                Yi[t | 6] = K = Yi[t | 6] ^ f, Yi[t | 7] = N = Yi[t | 7] ^ u, Yi[t | 8] = j = Yi[t | 8] ^ l, 
                Yi[t | 9] = D = Yi[t | 9] ^ c, Yi[t | 10] = G = Yi[t | 10] ^ w, Yi[t | 11] = F = Yi[t | 11] ^ y, 
                Yi[t | 12] = W = Yi[t | 12] ^ p, Yi[t | 13] = V = Yi[t | 13] ^ v, Yi[t | 14] = J = Yi[t | 14] ^ d, 
                Yi[t | 15] = Q = Yi[t | 15] ^ g;
                gn(_ ^ (z << 24 | O << 16 | T << 8 | R), k ^ (B << 24 | P << 16 | K << 8 | N), U ^ (j << 24 | D << 16 | G << 8 | F), x ^ (W << 24 | V << 16 | J << 8 | Q));
                Z = Z + 1 | 0;
                t = t + 16 | 0, e = e - 16 | 0, X = X + 16 | 0;
            }
            if ((e | 0) > 0) {
                en(r, b, m, S, A, E, L, H, I, q, M, C, Z >>> 24, Z >>> 16 & 255, Z >>> 8 & 255, Z & 255);
                z = Yi[t | 0] ^ i, O = (e | 0) > 1 ? Yi[t | 1] ^ n : 0, T = (e | 0) > 2 ? Yi[t | 2] ^ s : 0, 
                R = (e | 0) > 3 ? Yi[t | 3] ^ a : 0, B = (e | 0) > 4 ? Yi[t | 4] ^ h : 0, P = (e | 0) > 5 ? Yi[t | 5] ^ o : 0, 
                K = (e | 0) > 6 ? Yi[t | 6] ^ f : 0, N = (e | 0) > 7 ? Yi[t | 7] ^ u : 0, j = (e | 0) > 8 ? Yi[t | 8] ^ l : 0, 
                D = (e | 0) > 9 ? Yi[t | 9] ^ c : 0, G = (e | 0) > 10 ? Yi[t | 10] ^ w : 0, F = (e | 0) > 11 ? Yi[t | 11] ^ y : 0, 
                W = (e | 0) > 12 ? Yi[t | 12] ^ p : 0, V = (e | 0) > 13 ? Yi[t | 13] ^ v : 0, J = (e | 0) > 14 ? Yi[t | 14] ^ d : 0;
                Q = 0;
                Yi[t] = z;
                if ((e | 0) > 1) Yi[t | 1] = O;
                if ((e | 0) > 2) Yi[t | 2] = T;
                if ((e | 0) > 3) Yi[t | 3] = R;
                if ((e | 0) > 4) Yi[t | 4] = B;
                if ((e | 0) > 5) Yi[t | 5] = P;
                if ((e | 0) > 6) Yi[t | 6] = K;
                if ((e | 0) > 7) Yi[t | 7] = N;
                if ((e | 0) > 8) Yi[t | 8] = j;
                if ((e | 0) > 9) Yi[t | 9] = D;
                if ((e | 0) > 10) Yi[t | 10] = G;
                if ((e | 0) > 11) Yi[t | 11] = F;
                if ((e | 0) > 12) Yi[t | 12] = W;
                if ((e | 0) > 13) Yi[t | 13] = V;
                if ((e | 0) > 14) Yi[t | 14] = J;
                gn(_ ^ (z << 24 | O << 16 | T << 8 | R), k ^ (B << 24 | P << 16 | K << 8 | N), U ^ (j << 24 | D << 16 | G << 8 | F), x ^ (W << 24 | V << 16 | J << 8 | Q));
                Z = Z + 1 | 0;
                X = X + e | 0;
            }
            i = _ >>> 24, n = _ >>> 16 & 255, s = _ >>> 8 & 255, a = _ & 255, h = k >>> 24, 
            o = k >>> 16 & 255, f = k >>> 8 & 255, u = k & 255, l = U >>> 24, c = U >>> 16 & 255, 
            w = U >>> 8 & 255, y = U & 255, p = x >>> 24, v = x >>> 16 & 255, d = x >>> 8 & 255, 
            g = x & 255;
            return X | 0;
        }
        function An(t, e, r, b, m, S, A, E, L, H, I, q, M, C, Z) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            b = b | 0;
            m = m | 0;
            S = S | 0;
            A = A | 0;
            E = E | 0;
            L = L | 0;
            H = H | 0;
            I = I | 0;
            q = q | 0;
            M = M | 0;
            C = C | 0;
            Z = Z | 0;
            var z = 0, O = 0, T = 0, R = 0, B = 0, P = 0, K = 0, N = 0, j = 0, D = 0, G = 0, F = 0, W = 0, V = 0, J = 0, Q = 0, X = 0;
            if (t & 15) return -1;
            while ((e | 0) >= 16) {
                z = Yi[t | 0] | 0, O = Yi[t | 1] | 0, T = Yi[t | 2] | 0, R = Yi[t | 3] | 0, B = Yi[t | 4] | 0, 
                P = Yi[t | 5] | 0, K = Yi[t | 6] | 0, N = Yi[t | 7] | 0, j = Yi[t | 8] | 0, D = Yi[t | 9] | 0, 
                G = Yi[t | 10] | 0, F = Yi[t | 11] | 0, W = Yi[t | 12] | 0, V = Yi[t | 13] | 0, 
                J = Yi[t | 14] | 0, Q = Yi[t | 15] | 0;
                gn(_ ^ (z << 24 | O << 16 | T << 8 | R), k ^ (B << 24 | P << 16 | K << 8 | N), U ^ (j << 24 | D << 16 | G << 8 | F), x ^ (W << 24 | V << 16 | J << 8 | Q));
                en(r, b, m, S, A, E, L, H, I, q, M, C, Z >>> 24, Z >>> 16 & 255, Z >>> 8 & 255, Z & 255);
                Yi[t | 0] = z ^ i, Yi[t | 1] = O ^ n, Yi[t | 2] = T ^ s, Yi[t | 3] = R ^ a, Yi[t | 4] = B ^ h, 
                Yi[t | 5] = P ^ o, Yi[t | 6] = K ^ f, Yi[t | 7] = N ^ u, Yi[t | 8] = j ^ l, Yi[t | 9] = D ^ c, 
                Yi[t | 10] = G ^ w, Yi[t | 11] = F ^ y, Yi[t | 12] = W ^ p, Yi[t | 13] = V ^ v, 
                Yi[t | 14] = J ^ d, Yi[t | 15] = Q ^ g;
                Z = Z + 1 | 0;
                t = t + 16 | 0, e = e - 16 | 0, X = X + 16 | 0;
            }
            if ((e | 0) > 0) {
                z = Yi[t | 0] | 0, O = (e | 0) > 1 ? Yi[t | 1] | 0 : 0, T = (e | 0) > 2 ? Yi[t | 2] | 0 : 0, 
                R = (e | 0) > 3 ? Yi[t | 3] | 0 : 0, B = (e | 0) > 4 ? Yi[t | 4] | 0 : 0, P = (e | 0) > 5 ? Yi[t | 5] | 0 : 0, 
                K = (e | 0) > 6 ? Yi[t | 6] | 0 : 0, N = (e | 0) > 7 ? Yi[t | 7] | 0 : 0, j = (e | 0) > 8 ? Yi[t | 8] | 0 : 0, 
                D = (e | 0) > 9 ? Yi[t | 9] | 0 : 0, G = (e | 0) > 10 ? Yi[t | 10] | 0 : 0, F = (e | 0) > 11 ? Yi[t | 11] | 0 : 0, 
                W = (e | 0) > 12 ? Yi[t | 12] | 0 : 0, V = (e | 0) > 13 ? Yi[t | 13] | 0 : 0, J = (e | 0) > 14 ? Yi[t | 14] | 0 : 0;
                Q = 0;
                gn(_ ^ (z << 24 | O << 16 | T << 8 | R), k ^ (B << 24 | P << 16 | K << 8 | N), U ^ (j << 24 | D << 16 | G << 8 | F), x ^ (W << 24 | V << 16 | J << 8 | Q));
                en(r, b, m, S, A, E, L, H, I, q, M, C, Z >>> 24, Z >>> 16 & 255, Z >>> 8 & 255, Z & 255);
                Yi[t] = z ^ i;
                if ((e | 0) > 1) Yi[t | 1] = O ^ n;
                if ((e | 0) > 2) Yi[t | 2] = T ^ s;
                if ((e | 0) > 3) Yi[t | 3] = R ^ a;
                if ((e | 0) > 4) Yi[t | 4] = B ^ h;
                if ((e | 0) > 5) Yi[t | 5] = P ^ o;
                if ((e | 0) > 6) Yi[t | 6] = K ^ f;
                if ((e | 0) > 7) Yi[t | 7] = N ^ u;
                if ((e | 0) > 8) Yi[t | 8] = j ^ l;
                if ((e | 0) > 9) Yi[t | 9] = D ^ c;
                if ((e | 0) > 10) Yi[t | 10] = G ^ w;
                if ((e | 0) > 11) Yi[t | 11] = F ^ y;
                if ((e | 0) > 12) Yi[t | 12] = W ^ p;
                if ((e | 0) > 13) Yi[t | 13] = V ^ v;
                if ((e | 0) > 14) Yi[t | 14] = J ^ d;
                Z = Z + 1 | 0;
                X = X + e | 0;
            }
            i = _ >>> 24, n = _ >>> 16 & 255, s = _ >>> 8 & 255, a = _ & 255, h = k >>> 24, 
            o = k >>> 16 & 255, f = k >>> 8 & 255, u = k & 255, l = U >>> 24, c = U >>> 16 & 255, 
            w = U >>> 8 & 255, y = U & 255, p = x >>> 24, v = x >>> 16 & 255, d = x >>> 8 & 255, 
            g = x & 255;
            return X | 0;
        }
        return {
            init_state: nn,
            save_state: sn,
            init_key_128: an,
            init_key_256: hn,
            ecb_encrypt: on,
            ecb_decrypt: fn,
            cbc_encrypt: un,
            cbc_decrypt: ln,
            cbc_mac: cn,
            ctr_encrypt: wn,
            ctr_decrypt: wn,
            ccm_encrypt: yn,
            ccm_decrypt: pn,
            cfb_encrypt: vn,
            cfb_decrypt: dn,
            gcm_init: bn,
            gcm_ghash: mn,
            gcm_encrypt: Sn,
            gcm_decrypt: An
        };
    }
    function v(t, e, r) {
        return new Uint8Array(r).set(_t), p(t, e, r);
    }
    function d(t) {
        if ((t = t || {}).heapSize = t.heapSize || 4096, t.heapSize <= 0 || t.heapSize % 4096) throw new i("heapSize must be a positive number and multiple of 4096");
        this.BLOCK_SIZE = Ut, this.heap = t.heap || new Uint8Array(t.heapSize), this.asm = t.asm || v(e, null, this.heap.buffer), 
        this.pos = kt, this.len = 0, this.key = null, this.result = null, this.reset(t);
    }
    function g(t) {
        t = t || {}, this.result = null, this.pos = kt, this.len = 0;
        var e = this.asm, r = t.key;
        if (void 0 !== r) {
            if (c(r) || w(r)) r = new Uint8Array(r); else {
                if (!l(r)) throw new TypeError("unexpected key type");
                r = s(r);
            }
            if (16 === r.length) e.init_key_128.call(e, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11], r[12], r[13], r[14], r[15]); else {
                if (24 === r.length) throw new i("illegal key size");
                if (32 !== r.length) throw new i("illegal key size");
                e.init_key_256.call(e, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11], r[12], r[13], r[14], r[15], r[16], r[17], r[18], r[19], r[20], r[21], r[22], r[23], r[24], r[25], r[26], r[27], r[28], r[29], r[30], r[31]);
            }
            this.key = r;
        }
        return this;
    }
    function b(t, e, r, i, n) {
        var s = t.length - e, a = n > s ? s : n;
        return t.set(r.subarray(i, i + a), e), a;
    }
    function m(t) {
        this.padding = !0, this.mode = "cbc", this.iv = null, d.call(this, t);
    }
    function S(t) {
        t = t || {}, g.call(this, t);
        var e = t.padding;
        return this.padding = void 0 === e || !!e, function(t) {
            var e = this.asm;
            if (void 0 !== t) {
                if (c(t) || w(t)) t = new Uint8Array(t); else {
                    if (!l(t)) throw new TypeError("unexpected iv type");
                    t = s(t);
                }
                if (t.length !== Ut) throw new i("illegal iv size");
                this.iv = t, e.init_state.call(e, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10], t[11], t[12], t[13], t[14], t[15]);
            } else this.iv = null, e.init_state.call(e, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }.call(this, t.iv), this;
    }
    function A(t) {
        if (!this.key) throw new r("no key is associated with the instance");
        if (l(t) && (t = s(t)), c(t) && (t = new Uint8Array(t)), !w(t)) throw new TypeError("data isn't of expected type");
        for (var e = 0, i = t.length || 0, n = this.asm, a = this.heap, h = this.pos, o = this.len, f = 0, u = Ut * Math.floor((o + i) / Ut), y = 0, p = new Uint8Array(u); i > 0; ) o += y = b(a, h + o, t, e, i), 
        e += y, i -= y, y = n.cbc_encrypt(h, o), p.set(a.subarray(h, h + y), f), f += y, 
        o > y ? (h += y, o -= y) : (h = kt, o = 0);
        return this.result = p, this.pos = h, this.len = o, this;
    }
    function E() {
        if (!this.key) throw new r("no key is associated with the instance");
        var t = this.asm, e = this.heap, n = this.padding, s = this.pos, a = this.len, h = Ut * Math.ceil(a / Ut);
        if (a % Ut == 0) n && (h += Ut); else if (!n) throw new i("data length must be a multiple of " + Ut);
        var o = new Uint8Array(h);
        if (h > a) {
            for (var f = Ut - a % Ut, u = 0; f > u; ++u) e[s + a + u] = f;
            a += f;
        }
        return a > 0 && (t.cbc_encrypt(s, a), o.set(e.subarray(s, s + a))), this.result = o, 
        this.pos = kt, this.len = 0, this;
    }
    function _(t) {
        if (!this.key) throw new r("no key is associated with the instance");
        if (l(t) && (t = s(t)), c(t) && (t = new Uint8Array(t)), !w(t)) throw new TypeError("data isn't of expected type");
        for (var e = 0, i = t.length || 0, n = this.asm, a = this.heap, h = this.padding, o = this.pos, f = this.len, u = 0, y = Ut * Math.floor((f + i) / Ut), p = 0, v = new Uint8Array(y); i > 0; ) f += p = b(a, o + f, t, e, i), 
        e += p, i -= p, p = n.cbc_decrypt(o, f - (h && 0 === i && f % Ut == 0 ? Ut : 0)), 
        v.set(a.subarray(o, o + p), u), u += p, f > p ? (o += p, f -= p) : (o = kt, f = 0);
        return this.result = v.subarray(0, u), this.pos = o, this.len = f, this;
    }
    function k() {
        if (!this.key) throw new r("no key is associated with the instance");
        var t = this.asm, e = this.heap, n = this.padding, s = this.pos, a = this.len;
        if (0 === a) {
            if (n) throw new r("padding not found");
            return this.result = new Uint8Array(0), this.pos = kt, this.len = 0, this;
        }
        if (a % Ut != 0) throw new i("data length must be a multiple of " + Ut);
        var h = new Uint8Array(a);
        if (a > 0 && (t.cbc_decrypt(s, a), h.set(e.subarray(s, s + a))), n) {
            var o = h[a - 1];
            h = h.subarray(0, a - o);
        }
        return this.result = h, this.pos = kt, this.len = 0, this;
    }
    function U(t) {
        this.padding = !1, this.mode = "gcm", this.tagSize = Ut, this.adata = null, this.iv = null, 
        this.counter = 1, d.call(this, t);
    }
    function x(t) {
        for (var e = this.asm, r = this.heap, i = 0, n = t.length || 0, s = kt, a = 0, h = 0; n > 0; ) a += h = b(r, s + a, t, i, n), 
        i += h, n -= h, s += h = e.gcm_ghash(s, a), (a -= h) || (s = kt);
        if (a > 0) {
            for (;16 > a; ) r[s | a++] = 0;
            e.gcm_ghash(s, a);
        }
    }
    function L(t) {
        t = t || {};
        var e = this.asm, r = this.heap;
        g.call(this, t), e.gcm_init();
        var n = t.iv;
        if (void 0 !== n && null !== n) {
            if (c(n) || w(n)) n = new Uint8Array(n); else {
                if (!l(n)) throw new TypeError("unexpected iv type");
                n = s(n);
            }
            var a = n.length || 0;
            12 !== a ? (e.init_state(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), x.call(this, n), 
            r[0 | kt] = r[1 | kt] = r[2 | kt] = r[3 | kt] = r[4 | kt] = r[5 | kt] = r[6 | kt] = r[7 | kt] = r[8 | kt] = r[9 | kt] = r[10 | kt] = 0, 
            r[11 | kt] = a >>> 29, r[12 | kt] = a >>> 21 & 255, r[13 | kt] = a >>> 13 & 255, 
            r[14 | kt] = a >>> 5 & 255, r[15 | kt] = a << 3 & 255, e.gcm_ghash(kt, Ut), e.save_state(kt), 
            this.iv = new Uint8Array(r.subarray(kt, kt + Ut))) : (this.iv = new Uint8Array(16), 
            this.iv.set(n), this.iv[15] = 1);
        } else this.iv = new Uint8Array(16), this.iv[15] = 1;
        var h = t.counter;
        if (void 0 !== h) {
            if (!u(h)) throw new TypeError("counter must be a number");
            if (1 > h || h > 4294967295) throw new RangeError("counter must be a positive 32-bit integer");
            this.counter = h;
        } else this.counter = 1;
        var o = t.tagSize;
        if (void 0 !== o) {
            if (!u(o)) throw new TypeError("tagSize must be a number");
            if (4 > o || o > 16) throw new i("illegal tagSize value");
            this.tagSize = o;
        } else this.tagSize = Ut;
        e.init_state(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        var f = t.adata;
        if (void 0 !== f && null !== f) {
            if (c(f) || w(f)) f = new Uint8Array(f); else {
                if (!l(f)) throw new TypeError("unexpected adata type");
                f = s(f);
            }
            if (0 === f.length || f.length > It) throw new i("illegal adata length");
            x.call(this, f), this.adata = f;
        } else this.adata = null;
        return this;
    }
    function H(t) {
        if (!this.key) throw new r("no key is associated with the instance");
        if (l(t) && (t = s(t)), c(t) && (t = new Uint8Array(t)), !w(t)) throw new TypeError("data isn't of expected type");
        var e = this.asm, i = this.heap, n = this.iv, a = this.counter, h = 0, o = t.length || 0, f = this.pos, u = this.len, y = 0, p = Ut * Math.floor((u + o) / Ut), v = 0, d = new Uint8Array(p);
        if ((a - 1 << 4) + u + o > It) throw new r("counter overflow");
        for (;o > 0; ) {
            u += v = b(i, f + u, t, h, o), h += v, o -= v;
            var g = n[12] << 24 | n[13] << 16 | n[14] << 8 | n[15];
            (v = e.gcm_encrypt(f, -15 & u, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], g + a | 0)) && d.set(i.subarray(f, f + v), y), 
            a += v >>> 4, y += v, f += v, (u -= v) || (f = kt);
        }
        return this.result = d, this.counter = a, this.pos = f, this.len = u, this;
    }
    function I() {
        if (!this.key) throw new r("no key is associated with the instance");
        var t = this.asm, e = this.heap, i = this.iv, n = this.adata, s = this.counter, a = this.tagSize, h = this.pos, o = this.len, f = 0, u = new Uint8Array(o + a), l = i[12] << 24 | i[13] << 16 | i[14] << 8 | i[15];
        o > 0 && ((f = t.gcm_encrypt(h, o, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], l + s | 0)) && u.set(e.subarray(h, h + f)));
        var c = null !== n && n.length || 0, w = (s - 1 << 4) + f;
        return e[0 | kt] = e[1 | kt] = e[2 | kt] = 0, e[3 | kt] = c >>> 29, e[4 | kt] = c >>> 21, 
        e[5 | kt] = c >>> 13 & 255, e[6 | kt] = c >>> 5 & 255, e[7 | kt] = c << 3 & 255, 
        e[8 | kt] = e[9 | kt] = e[10 | kt] = 0, e[11 | kt] = w >>> 29, e[12 | kt] = w >>> 21 & 255, 
        e[13 | kt] = w >>> 13 & 255, e[14 | kt] = w >>> 5 & 255, e[15 | kt] = w << 3 & 255, 
        t.gcm_ghash(kt, Ut), t.save_state(kt), t.gcm_encrypt(kt, Ut, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], l), 
        u.set(e.subarray(kt, kt + a), f), this.result = u, this.counter = 1, this.pos = kt, 
        this.len = 0, this;
    }
    function q(t) {
        if (!this.key) throw new r("no key is associated with the instance");
        if (l(t) && (t = s(t)), c(t) && (t = new Uint8Array(t)), !w(t)) throw new TypeError("data isn't of expected type");
        var e = this.asm, i = this.heap, n = this.iv, a = this.counter, h = this.tagSize, o = 0, f = t.length || 0, u = this.pos, y = this.len, p = 0, v = Ut * Math.floor((y + f - h) / Ut), d = 0, g = new Uint8Array(v);
        if ((a - 1 << 4) + y + f - h > It) throw new r("counter overflow");
        for (;f > 0; ) {
            y += d = b(i, u + y, t, o, f), o += d, f -= d;
            var m = n[12] << 24 | n[13] << 16 | n[14] << 8 | n[15];
            (d = e.gcm_decrypt(u, -15 & Math.min(y, y + f - h), n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], m + a | 0)) && g.set(i.subarray(u, u + d), p), 
            a += d >>> 4, p += d, u += d, (y -= d) || (u = kt);
        }
        return this.result = g, this.counter = a, this.pos = u, this.len = y, this;
    }
    function M() {
        if (!this.key) throw new r("no key is associated with the instance");
        var t = this.asm, e = this.heap, i = this.iv, s = this.adata, a = this.counter, h = this.tagSize, o = this.pos, f = this.len, u = f - h, l = 0, c = new Uint8Array(u), w = new Uint8Array(e.subarray(o + u, o + f)), y = i[12] << 24 | i[13] << 16 | i[14] << 8 | i[15];
        f > 0 && ((l = t.gcm_decrypt(o, u, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], y + a | 0)) && c.set(e.subarray(o, o + l)));
        var p = null !== s && s.length || 0, v = (a - 1 << 4) + l;
        e[0 | kt] = e[1 | kt] = e[2 | kt] = 0, e[3 | kt] = p >>> 29, e[4 | kt] = p >>> 21, 
        e[5 | kt] = p >>> 13 & 255, e[6 | kt] = p >>> 5 & 255, e[7 | kt] = p << 3 & 255, 
        e[8 | kt] = e[9 | kt] = e[10 | kt] = 0, e[11 | kt] = v >>> 29, e[12 | kt] = v >>> 21 & 255, 
        e[13 | kt] = v >>> 13 & 255, e[14 | kt] = v >>> 5 & 255, e[15 | kt] = v << 3 & 255, 
        t.gcm_ghash(kt, Ut), t.save_state(kt), t.gcm_encrypt(kt, Ut, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], y);
        for (var d = 0, g = 0; h > g; ++g) d |= w[g] ^ e[kt | g];
        if (d) throw new n("data integrity check failed");
        return this.result = c, this.counter = 1, this.pos = kt, this.len = 0, this;
    }
    function C(t) {
        if ((t = t || {}).heapSize = t.heapSize || 4096, t.heapSize <= 0 || t.heapSize % 4096) throw new i("heapSize must be a positive number and multiple of 4096");
        this.heap = t.heap || new Uint8Array(t.heapSize), this.asm = t.asm || function(t, e, r) {
            "use asm";
            var i = 0, n = 0, s = 0, a = 0, h = 0, o = 0, f = 0, u = 0, l = 0, c = 0, w = 0, y = 0, p = 0, v = 0, d = 0, g = 0, b = new t.Uint8Array(r);
            function m(t, e, r, o, f, u, l, c, w, y, p, v, d, g, b, m) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                o = o | 0;
                f = f | 0;
                u = u | 0;
                l = l | 0;
                c = c | 0;
                w = w | 0;
                y = y | 0;
                p = p | 0;
                v = v | 0;
                d = d | 0;
                g = g | 0;
                b = b | 0;
                m = m | 0;
                var S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0, L = 0, H = 0, I = 0, q = 0, M = 0, C = 0, Z = 0, z = 0, O = 0, T = 0, R = 0, B = 0, P = 0, K = 0, N = 0, j = 0, D = 0, G = 0, F = 0, W = 0, V = 0, J = 0, Q = 0, X = 0, Y = 0, $ = 0, tt = 0, et = 0, rt = 0, it = 0, nt = 0, st = 0, at = 0, ht = 0, ot = 0, ft = 0, ut = 0, lt = 0, ct = 0, wt = 0, yt = 0, pt = 0, vt = 0, dt = 0, gt = 0, bt = 0, mt = 0, St = 0, At = 0, Et = 0, _t = 0, kt = 0, Ut = 0, xt = 0, Lt = 0, Ht = 0, It = 0, qt = 0, Mt = 0, Ct = 0, Zt = 0, zt = 0, Ot = 0, Tt = 0;
                S = i;
                A = n;
                E = s;
                _ = a;
                k = h;
                x = t + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = e + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = r + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = o + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = f + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = u + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = l + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = c + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = w + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = y + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = p + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = v + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = d + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = g + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = b + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                x = m + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = g ^ w ^ r ^ t;
                L = U << 1 | U >>> 31;
                x = L + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = b ^ y ^ o ^ e;
                H = U << 1 | U >>> 31;
                x = H + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = m ^ p ^ f ^ r;
                I = U << 1 | U >>> 31;
                x = I + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = L ^ v ^ u ^ o;
                q = U << 1 | U >>> 31;
                x = q + (S << 5 | S >>> 27) + k + (A & E | ~A & _) + 1518500249 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = H ^ d ^ l ^ f;
                M = U << 1 | U >>> 31;
                x = M + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = I ^ g ^ c ^ u;
                C = U << 1 | U >>> 31;
                x = C + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = q ^ b ^ w ^ l;
                Z = U << 1 | U >>> 31;
                x = Z + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = M ^ m ^ y ^ c;
                z = U << 1 | U >>> 31;
                x = z + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = C ^ L ^ p ^ w;
                O = U << 1 | U >>> 31;
                x = O + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Z ^ H ^ v ^ y;
                T = U << 1 | U >>> 31;
                x = T + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = z ^ I ^ d ^ p;
                R = U << 1 | U >>> 31;
                x = R + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = O ^ q ^ g ^ v;
                B = U << 1 | U >>> 31;
                x = B + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = T ^ M ^ b ^ d;
                P = U << 1 | U >>> 31;
                x = P + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = R ^ C ^ m ^ g;
                K = U << 1 | U >>> 31;
                x = K + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = B ^ Z ^ L ^ b;
                N = U << 1 | U >>> 31;
                x = N + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = P ^ z ^ H ^ m;
                j = U << 1 | U >>> 31;
                x = j + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = K ^ O ^ I ^ L;
                D = U << 1 | U >>> 31;
                x = D + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = N ^ T ^ q ^ H;
                G = U << 1 | U >>> 31;
                x = G + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = j ^ R ^ M ^ I;
                F = U << 1 | U >>> 31;
                x = F + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = D ^ B ^ C ^ q;
                W = U << 1 | U >>> 31;
                x = W + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = G ^ P ^ Z ^ M;
                V = U << 1 | U >>> 31;
                x = V + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = F ^ K ^ z ^ C;
                J = U << 1 | U >>> 31;
                x = J + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = W ^ N ^ O ^ Z;
                Q = U << 1 | U >>> 31;
                x = Q + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = V ^ j ^ T ^ z;
                X = U << 1 | U >>> 31;
                x = X + (S << 5 | S >>> 27) + k + (A ^ E ^ _) + 1859775393 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = J ^ D ^ R ^ O;
                Y = U << 1 | U >>> 31;
                x = Y + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Q ^ G ^ B ^ T;
                $ = U << 1 | U >>> 31;
                x = $ + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = X ^ F ^ P ^ R;
                tt = U << 1 | U >>> 31;
                x = tt + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Y ^ W ^ K ^ B;
                et = U << 1 | U >>> 31;
                x = et + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = $ ^ V ^ N ^ P;
                rt = U << 1 | U >>> 31;
                x = rt + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = tt ^ J ^ j ^ K;
                it = U << 1 | U >>> 31;
                x = it + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = et ^ Q ^ D ^ N;
                nt = U << 1 | U >>> 31;
                x = nt + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = rt ^ X ^ G ^ j;
                st = U << 1 | U >>> 31;
                x = st + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = it ^ Y ^ F ^ D;
                at = U << 1 | U >>> 31;
                x = at + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = nt ^ $ ^ W ^ G;
                ht = U << 1 | U >>> 31;
                x = ht + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = st ^ tt ^ V ^ F;
                ot = U << 1 | U >>> 31;
                x = ot + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = at ^ et ^ J ^ W;
                ft = U << 1 | U >>> 31;
                x = ft + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = ht ^ rt ^ Q ^ V;
                ut = U << 1 | U >>> 31;
                x = ut + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = ot ^ it ^ X ^ J;
                lt = U << 1 | U >>> 31;
                x = lt + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = ft ^ nt ^ Y ^ Q;
                ct = U << 1 | U >>> 31;
                x = ct + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = ut ^ st ^ $ ^ X;
                wt = U << 1 | U >>> 31;
                x = wt + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = lt ^ at ^ tt ^ Y;
                yt = U << 1 | U >>> 31;
                x = yt + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = ct ^ ht ^ et ^ $;
                pt = U << 1 | U >>> 31;
                x = pt + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = wt ^ ot ^ rt ^ tt;
                vt = U << 1 | U >>> 31;
                x = vt + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = yt ^ ft ^ it ^ et;
                dt = U << 1 | U >>> 31;
                x = dt + (S << 5 | S >>> 27) + k + (A & E | A & _ | E & _) - 1894007588 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = pt ^ ut ^ nt ^ rt;
                gt = U << 1 | U >>> 31;
                x = gt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = vt ^ lt ^ st ^ it;
                bt = U << 1 | U >>> 31;
                x = bt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = dt ^ ct ^ at ^ nt;
                mt = U << 1 | U >>> 31;
                x = mt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = gt ^ wt ^ ht ^ st;
                St = U << 1 | U >>> 31;
                x = St + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = bt ^ yt ^ ot ^ at;
                At = U << 1 | U >>> 31;
                x = At + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = mt ^ pt ^ ft ^ ht;
                Et = U << 1 | U >>> 31;
                x = Et + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = St ^ vt ^ ut ^ ot;
                _t = U << 1 | U >>> 31;
                x = _t + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = At ^ dt ^ lt ^ ft;
                kt = U << 1 | U >>> 31;
                x = kt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Et ^ gt ^ ct ^ ut;
                Ut = U << 1 | U >>> 31;
                x = Ut + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = _t ^ bt ^ wt ^ lt;
                xt = U << 1 | U >>> 31;
                x = xt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = kt ^ mt ^ yt ^ ct;
                Lt = U << 1 | U >>> 31;
                x = Lt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Ut ^ St ^ pt ^ wt;
                Ht = U << 1 | U >>> 31;
                x = Ht + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = xt ^ At ^ vt ^ yt;
                It = U << 1 | U >>> 31;
                x = It + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Lt ^ Et ^ dt ^ pt;
                qt = U << 1 | U >>> 31;
                x = qt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Ht ^ _t ^ gt ^ vt;
                Mt = U << 1 | U >>> 31;
                x = Mt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = It ^ kt ^ bt ^ dt;
                Ct = U << 1 | U >>> 31;
                x = Ct + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = qt ^ Ut ^ mt ^ gt;
                Zt = U << 1 | U >>> 31;
                x = Zt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Mt ^ xt ^ St ^ bt;
                zt = U << 1 | U >>> 31;
                x = zt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Ct ^ Lt ^ At ^ mt;
                Ot = U << 1 | U >>> 31;
                x = Ot + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                U = Zt ^ Ht ^ Et ^ St;
                Tt = U << 1 | U >>> 31;
                x = Tt + (S << 5 | S >>> 27) + k + (A ^ E ^ _) - 899497514 | 0;
                k = _;
                _ = E;
                E = A << 30 | A >>> 2;
                A = S;
                S = x;
                i = i + S | 0;
                n = n + A | 0;
                s = s + E | 0;
                a = a + _ | 0;
                h = h + k | 0;
            }
            function S(t) {
                t = t | 0;
                m(b[t | 0] << 24 | b[t | 1] << 16 | b[t | 2] << 8 | b[t | 3], b[t | 4] << 24 | b[t | 5] << 16 | b[t | 6] << 8 | b[t | 7], b[t | 8] << 24 | b[t | 9] << 16 | b[t | 10] << 8 | b[t | 11], b[t | 12] << 24 | b[t | 13] << 16 | b[t | 14] << 8 | b[t | 15], b[t | 16] << 24 | b[t | 17] << 16 | b[t | 18] << 8 | b[t | 19], b[t | 20] << 24 | b[t | 21] << 16 | b[t | 22] << 8 | b[t | 23], b[t | 24] << 24 | b[t | 25] << 16 | b[t | 26] << 8 | b[t | 27], b[t | 28] << 24 | b[t | 29] << 16 | b[t | 30] << 8 | b[t | 31], b[t | 32] << 24 | b[t | 33] << 16 | b[t | 34] << 8 | b[t | 35], b[t | 36] << 24 | b[t | 37] << 16 | b[t | 38] << 8 | b[t | 39], b[t | 40] << 24 | b[t | 41] << 16 | b[t | 42] << 8 | b[t | 43], b[t | 44] << 24 | b[t | 45] << 16 | b[t | 46] << 8 | b[t | 47], b[t | 48] << 24 | b[t | 49] << 16 | b[t | 50] << 8 | b[t | 51], b[t | 52] << 24 | b[t | 53] << 16 | b[t | 54] << 8 | b[t | 55], b[t | 56] << 24 | b[t | 57] << 16 | b[t | 58] << 8 | b[t | 59], b[t | 60] << 24 | b[t | 61] << 16 | b[t | 62] << 8 | b[t | 63]);
            }
            function A(t) {
                t = t | 0;
                b[t | 0] = i >>> 24;
                b[t | 1] = i >>> 16 & 255;
                b[t | 2] = i >>> 8 & 255;
                b[t | 3] = i & 255;
                b[t | 4] = n >>> 24;
                b[t | 5] = n >>> 16 & 255;
                b[t | 6] = n >>> 8 & 255;
                b[t | 7] = n & 255;
                b[t | 8] = s >>> 24;
                b[t | 9] = s >>> 16 & 255;
                b[t | 10] = s >>> 8 & 255;
                b[t | 11] = s & 255;
                b[t | 12] = a >>> 24;
                b[t | 13] = a >>> 16 & 255;
                b[t | 14] = a >>> 8 & 255;
                b[t | 15] = a & 255;
                b[t | 16] = h >>> 24;
                b[t | 17] = h >>> 16 & 255;
                b[t | 18] = h >>> 8 & 255;
                b[t | 19] = h & 255;
            }
            function E() {
                i = 1732584193;
                n = 4023233417;
                s = 2562383102;
                a = 271733878;
                h = 3285377520;
                o = 0;
            }
            function _(t, e, r, f, u, l) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                f = f | 0;
                u = u | 0;
                l = l | 0;
                i = t;
                n = e;
                s = r;
                a = f;
                h = u;
                o = l;
            }
            function k(t, e) {
                t = t | 0;
                e = e | 0;
                var r = 0;
                if (t & 63) return -1;
                while ((e | 0) >= 64) {
                    S(t);
                    t = t + 64 | 0;
                    e = e - 64 | 0;
                    r = r + 64 | 0;
                }
                o = o + r | 0;
                return r | 0;
            }
            function U(t, e, r) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                var i = 0, n = 0;
                if (t & 63) return -1;
                if (~r) if (r & 31) return -1;
                if ((e | 0) >= 64) {
                    i = k(t, e) | 0;
                    if ((i | 0) == -1) return -1;
                    t = t + i | 0;
                    e = e - i | 0;
                }
                i = i + e | 0;
                o = o + e | 0;
                b[t | e] = 128;
                if ((e | 0) >= 56) {
                    for (n = e + 1 | 0; (n | 0) < 64; n = n + 1 | 0) b[t | n] = 0;
                    S(t);
                    e = 0;
                    b[t | 0] = 0;
                }
                for (n = e + 1 | 0; (n | 0) < 59; n = n + 1 | 0) b[t | n] = 0;
                b[t | 59] = o >>> 29;
                b[t | 60] = o >>> 21 & 255;
                b[t | 61] = o >>> 13 & 255;
                b[t | 62] = o >>> 5 & 255;
                b[t | 63] = o << 3 & 255;
                S(t);
                if (~r) A(r);
                return i | 0;
            }
            function x() {
                i = f;
                n = u;
                s = l;
                a = c;
                h = w;
                o = 64;
            }
            function L() {
                i = y;
                n = p;
                s = v;
                a = d;
                h = g;
                o = 64;
            }
            function H(t, e, r, b, S, A, _, k, U, x, L, H, I, q, M, C) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                b = b | 0;
                S = S | 0;
                A = A | 0;
                _ = _ | 0;
                k = k | 0;
                U = U | 0;
                x = x | 0;
                L = L | 0;
                H = H | 0;
                I = I | 0;
                q = q | 0;
                M = M | 0;
                C = C | 0;
                E();
                m(t ^ 1549556828, e ^ 1549556828, r ^ 1549556828, b ^ 1549556828, S ^ 1549556828, A ^ 1549556828, _ ^ 1549556828, k ^ 1549556828, U ^ 1549556828, x ^ 1549556828, L ^ 1549556828, H ^ 1549556828, I ^ 1549556828, q ^ 1549556828, M ^ 1549556828, C ^ 1549556828);
                y = i;
                p = n;
                v = s;
                d = a;
                g = h;
                E();
                m(t ^ 909522486, e ^ 909522486, r ^ 909522486, b ^ 909522486, S ^ 909522486, A ^ 909522486, _ ^ 909522486, k ^ 909522486, U ^ 909522486, x ^ 909522486, L ^ 909522486, H ^ 909522486, I ^ 909522486, q ^ 909522486, M ^ 909522486, C ^ 909522486);
                f = i;
                u = n;
                l = s;
                c = a;
                w = h;
                o = 64;
            }
            function I(t, e, r) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                var o = 0, f = 0, u = 0, l = 0, c = 0, w = 0;
                if (t & 63) return -1;
                if (~r) if (r & 31) return -1;
                w = U(t, e, -1) | 0;
                o = i, f = n, u = s, l = a, c = h;
                L();
                m(o, f, u, l, c, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 672);
                if (~r) A(r);
                return w | 0;
            }
            function q(t, e, r, o, f) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                o = o | 0;
                f = f | 0;
                var u = 0, l = 0, c = 0, w = 0, y = 0, p = 0, v = 0, d = 0, g = 0, S = 0;
                if (t & 63) return -1;
                if (~f) if (f & 31) return -1;
                b[t + e | 0] = r >>> 24;
                b[t + e + 1 | 0] = r >>> 16 & 255;
                b[t + e + 2 | 0] = r >>> 8 & 255;
                b[t + e + 3 | 0] = r & 255;
                I(t, e + 4 | 0, -1) | 0;
                u = p = i, l = v = n, c = d = s, w = g = a, y = S = h;
                o = o - 1 | 0;
                while ((o | 0) > 0) {
                    x();
                    m(p, v, d, g, S, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 672);
                    p = i, v = n, d = s, g = a, S = h;
                    L();
                    m(p, v, d, g, S, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 672);
                    p = i, v = n, d = s, g = a, S = h;
                    u = u ^ i;
                    l = l ^ n;
                    c = c ^ s;
                    w = w ^ a;
                    y = y ^ h;
                    o = o - 1 | 0;
                }
                i = u;
                n = l;
                s = c;
                a = w;
                h = y;
                if (~f) A(f);
                return 0;
            }
            return {
                reset: E,
                init: _,
                process: k,
                finish: U,
                hmac_reset: x,
                hmac_init: H,
                hmac_finish: I,
                pbkdf2_generate_block: q
            };
        }(e, null, this.heap.buffer), this.BLOCK_SIZE = Tt, this.HASH_SIZE = Rt, this.reset();
    }
    function Z() {
        return null === Pt && (Pt = new C({
            heapSize: 1048576
        })), Pt;
    }
    function z(t) {
        if (void 0 === t) throw new SyntaxError("data required");
        return Z().reset().process(t).finish().result;
    }
    function O(t) {
        if ((t = t || {}).heapSize = t.heapSize || 4096, t.heapSize <= 0 || t.heapSize % 4096) throw new i("heapSize must be a positive number and multiple of 4096");
        this.heap = t.heap || new Uint8Array(t.heapSize), this.asm = t.asm || function(t, e, r) {
            "use asm";
            var i = 0, n = 0, s = 0, a = 0, h = 0, o = 0, f = 0, u = 0, l = 0, c = 0, w = 0, y = 0, p = 0, v = 0, d = 0, g = 0, b = 0, m = 0, S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0, L = new t.Uint8Array(r);
            function H(t, e, r, l, c, w, y, p, v, d, g, b, m, S, A, E) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                l = l | 0;
                c = c | 0;
                w = w | 0;
                y = y | 0;
                p = p | 0;
                v = v | 0;
                d = d | 0;
                g = g | 0;
                b = b | 0;
                m = m | 0;
                S = S | 0;
                A = A | 0;
                E = E | 0;
                var _ = 0, k = 0, U = 0, x = 0, L = 0, H = 0, I = 0, q = 0, M = 0;
                _ = i;
                k = n;
                U = s;
                x = a;
                L = h;
                H = o;
                I = f;
                q = u;
                M = t + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1116352408 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = e + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1899447441 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = r + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3049323471 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = l + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3921009573 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = c + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 961987163 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = w + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1508970993 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = y + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2453635748 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = p + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2870763221 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = v + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3624381080 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = d + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 310598401 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = g + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 607225278 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = b + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1426881987 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = m + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1925078388 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = S + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2162078206 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = A + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2614888103 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                M = E + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3248222580 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                t = M = (e >>> 7 ^ e >>> 18 ^ e >>> 3 ^ e << 25 ^ e << 14) + (A >>> 17 ^ A >>> 19 ^ A >>> 10 ^ A << 15 ^ A << 13) + t + d | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3835390401 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                e = M = (r >>> 7 ^ r >>> 18 ^ r >>> 3 ^ r << 25 ^ r << 14) + (E >>> 17 ^ E >>> 19 ^ E >>> 10 ^ E << 15 ^ E << 13) + e + g | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 4022224774 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                r = M = (l >>> 7 ^ l >>> 18 ^ l >>> 3 ^ l << 25 ^ l << 14) + (t >>> 17 ^ t >>> 19 ^ t >>> 10 ^ t << 15 ^ t << 13) + r + b | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 264347078 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                l = M = (c >>> 7 ^ c >>> 18 ^ c >>> 3 ^ c << 25 ^ c << 14) + (e >>> 17 ^ e >>> 19 ^ e >>> 10 ^ e << 15 ^ e << 13) + l + m | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 604807628 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                c = M = (w >>> 7 ^ w >>> 18 ^ w >>> 3 ^ w << 25 ^ w << 14) + (r >>> 17 ^ r >>> 19 ^ r >>> 10 ^ r << 15 ^ r << 13) + c + S | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 770255983 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                w = M = (y >>> 7 ^ y >>> 18 ^ y >>> 3 ^ y << 25 ^ y << 14) + (l >>> 17 ^ l >>> 19 ^ l >>> 10 ^ l << 15 ^ l << 13) + w + A | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1249150122 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                y = M = (p >>> 7 ^ p >>> 18 ^ p >>> 3 ^ p << 25 ^ p << 14) + (c >>> 17 ^ c >>> 19 ^ c >>> 10 ^ c << 15 ^ c << 13) + y + E | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1555081692 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                p = M = (v >>> 7 ^ v >>> 18 ^ v >>> 3 ^ v << 25 ^ v << 14) + (w >>> 17 ^ w >>> 19 ^ w >>> 10 ^ w << 15 ^ w << 13) + p + t | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1996064986 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                v = M = (d >>> 7 ^ d >>> 18 ^ d >>> 3 ^ d << 25 ^ d << 14) + (y >>> 17 ^ y >>> 19 ^ y >>> 10 ^ y << 15 ^ y << 13) + v + e | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2554220882 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                d = M = (g >>> 7 ^ g >>> 18 ^ g >>> 3 ^ g << 25 ^ g << 14) + (p >>> 17 ^ p >>> 19 ^ p >>> 10 ^ p << 15 ^ p << 13) + d + r | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2821834349 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                g = M = (b >>> 7 ^ b >>> 18 ^ b >>> 3 ^ b << 25 ^ b << 14) + (v >>> 17 ^ v >>> 19 ^ v >>> 10 ^ v << 15 ^ v << 13) + g + l | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2952996808 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                b = M = (m >>> 7 ^ m >>> 18 ^ m >>> 3 ^ m << 25 ^ m << 14) + (d >>> 17 ^ d >>> 19 ^ d >>> 10 ^ d << 15 ^ d << 13) + b + c | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3210313671 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                m = M = (S >>> 7 ^ S >>> 18 ^ S >>> 3 ^ S << 25 ^ S << 14) + (g >>> 17 ^ g >>> 19 ^ g >>> 10 ^ g << 15 ^ g << 13) + m + w | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3336571891 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                S = M = (A >>> 7 ^ A >>> 18 ^ A >>> 3 ^ A << 25 ^ A << 14) + (b >>> 17 ^ b >>> 19 ^ b >>> 10 ^ b << 15 ^ b << 13) + S + y | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3584528711 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                A = M = (E >>> 7 ^ E >>> 18 ^ E >>> 3 ^ E << 25 ^ E << 14) + (m >>> 17 ^ m >>> 19 ^ m >>> 10 ^ m << 15 ^ m << 13) + A + p | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 113926993 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                E = M = (t >>> 7 ^ t >>> 18 ^ t >>> 3 ^ t << 25 ^ t << 14) + (S >>> 17 ^ S >>> 19 ^ S >>> 10 ^ S << 15 ^ S << 13) + E + v | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 338241895 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                t = M = (e >>> 7 ^ e >>> 18 ^ e >>> 3 ^ e << 25 ^ e << 14) + (A >>> 17 ^ A >>> 19 ^ A >>> 10 ^ A << 15 ^ A << 13) + t + d | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 666307205 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                e = M = (r >>> 7 ^ r >>> 18 ^ r >>> 3 ^ r << 25 ^ r << 14) + (E >>> 17 ^ E >>> 19 ^ E >>> 10 ^ E << 15 ^ E << 13) + e + g | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 773529912 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                r = M = (l >>> 7 ^ l >>> 18 ^ l >>> 3 ^ l << 25 ^ l << 14) + (t >>> 17 ^ t >>> 19 ^ t >>> 10 ^ t << 15 ^ t << 13) + r + b | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1294757372 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                l = M = (c >>> 7 ^ c >>> 18 ^ c >>> 3 ^ c << 25 ^ c << 14) + (e >>> 17 ^ e >>> 19 ^ e >>> 10 ^ e << 15 ^ e << 13) + l + m | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1396182291 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                c = M = (w >>> 7 ^ w >>> 18 ^ w >>> 3 ^ w << 25 ^ w << 14) + (r >>> 17 ^ r >>> 19 ^ r >>> 10 ^ r << 15 ^ r << 13) + c + S | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1695183700 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                w = M = (y >>> 7 ^ y >>> 18 ^ y >>> 3 ^ y << 25 ^ y << 14) + (l >>> 17 ^ l >>> 19 ^ l >>> 10 ^ l << 15 ^ l << 13) + w + A | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1986661051 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                y = M = (p >>> 7 ^ p >>> 18 ^ p >>> 3 ^ p << 25 ^ p << 14) + (c >>> 17 ^ c >>> 19 ^ c >>> 10 ^ c << 15 ^ c << 13) + y + E | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2177026350 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                p = M = (v >>> 7 ^ v >>> 18 ^ v >>> 3 ^ v << 25 ^ v << 14) + (w >>> 17 ^ w >>> 19 ^ w >>> 10 ^ w << 15 ^ w << 13) + p + t | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2456956037 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                v = M = (d >>> 7 ^ d >>> 18 ^ d >>> 3 ^ d << 25 ^ d << 14) + (y >>> 17 ^ y >>> 19 ^ y >>> 10 ^ y << 15 ^ y << 13) + v + e | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2730485921 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                d = M = (g >>> 7 ^ g >>> 18 ^ g >>> 3 ^ g << 25 ^ g << 14) + (p >>> 17 ^ p >>> 19 ^ p >>> 10 ^ p << 15 ^ p << 13) + d + r | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2820302411 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                g = M = (b >>> 7 ^ b >>> 18 ^ b >>> 3 ^ b << 25 ^ b << 14) + (v >>> 17 ^ v >>> 19 ^ v >>> 10 ^ v << 15 ^ v << 13) + g + l | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3259730800 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                b = M = (m >>> 7 ^ m >>> 18 ^ m >>> 3 ^ m << 25 ^ m << 14) + (d >>> 17 ^ d >>> 19 ^ d >>> 10 ^ d << 15 ^ d << 13) + b + c | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3345764771 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                m = M = (S >>> 7 ^ S >>> 18 ^ S >>> 3 ^ S << 25 ^ S << 14) + (g >>> 17 ^ g >>> 19 ^ g >>> 10 ^ g << 15 ^ g << 13) + m + w | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3516065817 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                S = M = (A >>> 7 ^ A >>> 18 ^ A >>> 3 ^ A << 25 ^ A << 14) + (b >>> 17 ^ b >>> 19 ^ b >>> 10 ^ b << 15 ^ b << 13) + S + y | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3600352804 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                A = M = (E >>> 7 ^ E >>> 18 ^ E >>> 3 ^ E << 25 ^ E << 14) + (m >>> 17 ^ m >>> 19 ^ m >>> 10 ^ m << 15 ^ m << 13) + A + p | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 4094571909 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                E = M = (t >>> 7 ^ t >>> 18 ^ t >>> 3 ^ t << 25 ^ t << 14) + (S >>> 17 ^ S >>> 19 ^ S >>> 10 ^ S << 15 ^ S << 13) + E + v | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 275423344 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                t = M = (e >>> 7 ^ e >>> 18 ^ e >>> 3 ^ e << 25 ^ e << 14) + (A >>> 17 ^ A >>> 19 ^ A >>> 10 ^ A << 15 ^ A << 13) + t + d | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 430227734 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                e = M = (r >>> 7 ^ r >>> 18 ^ r >>> 3 ^ r << 25 ^ r << 14) + (E >>> 17 ^ E >>> 19 ^ E >>> 10 ^ E << 15 ^ E << 13) + e + g | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 506948616 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                r = M = (l >>> 7 ^ l >>> 18 ^ l >>> 3 ^ l << 25 ^ l << 14) + (t >>> 17 ^ t >>> 19 ^ t >>> 10 ^ t << 15 ^ t << 13) + r + b | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 659060556 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                l = M = (c >>> 7 ^ c >>> 18 ^ c >>> 3 ^ c << 25 ^ c << 14) + (e >>> 17 ^ e >>> 19 ^ e >>> 10 ^ e << 15 ^ e << 13) + l + m | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 883997877 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                c = M = (w >>> 7 ^ w >>> 18 ^ w >>> 3 ^ w << 25 ^ w << 14) + (r >>> 17 ^ r >>> 19 ^ r >>> 10 ^ r << 15 ^ r << 13) + c + S | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 958139571 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                w = M = (y >>> 7 ^ y >>> 18 ^ y >>> 3 ^ y << 25 ^ y << 14) + (l >>> 17 ^ l >>> 19 ^ l >>> 10 ^ l << 15 ^ l << 13) + w + A | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1322822218 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                y = M = (p >>> 7 ^ p >>> 18 ^ p >>> 3 ^ p << 25 ^ p << 14) + (c >>> 17 ^ c >>> 19 ^ c >>> 10 ^ c << 15 ^ c << 13) + y + E | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1537002063 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                p = M = (v >>> 7 ^ v >>> 18 ^ v >>> 3 ^ v << 25 ^ v << 14) + (w >>> 17 ^ w >>> 19 ^ w >>> 10 ^ w << 15 ^ w << 13) + p + t | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1747873779 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                v = M = (d >>> 7 ^ d >>> 18 ^ d >>> 3 ^ d << 25 ^ d << 14) + (y >>> 17 ^ y >>> 19 ^ y >>> 10 ^ y << 15 ^ y << 13) + v + e | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 1955562222 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                d = M = (g >>> 7 ^ g >>> 18 ^ g >>> 3 ^ g << 25 ^ g << 14) + (p >>> 17 ^ p >>> 19 ^ p >>> 10 ^ p << 15 ^ p << 13) + d + r | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2024104815 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                g = M = (b >>> 7 ^ b >>> 18 ^ b >>> 3 ^ b << 25 ^ b << 14) + (v >>> 17 ^ v >>> 19 ^ v >>> 10 ^ v << 15 ^ v << 13) + g + l | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2227730452 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                b = M = (m >>> 7 ^ m >>> 18 ^ m >>> 3 ^ m << 25 ^ m << 14) + (d >>> 17 ^ d >>> 19 ^ d >>> 10 ^ d << 15 ^ d << 13) + b + c | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2361852424 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                m = M = (S >>> 7 ^ S >>> 18 ^ S >>> 3 ^ S << 25 ^ S << 14) + (g >>> 17 ^ g >>> 19 ^ g >>> 10 ^ g << 15 ^ g << 13) + m + w | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2428436474 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                S = M = (A >>> 7 ^ A >>> 18 ^ A >>> 3 ^ A << 25 ^ A << 14) + (b >>> 17 ^ b >>> 19 ^ b >>> 10 ^ b << 15 ^ b << 13) + S + y | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 2756734187 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                A = M = (E >>> 7 ^ E >>> 18 ^ E >>> 3 ^ E << 25 ^ E << 14) + (m >>> 17 ^ m >>> 19 ^ m >>> 10 ^ m << 15 ^ m << 13) + A + p | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3204031479 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                E = M = (t >>> 7 ^ t >>> 18 ^ t >>> 3 ^ t << 25 ^ t << 14) + (S >>> 17 ^ S >>> 19 ^ S >>> 10 ^ S << 15 ^ S << 13) + E + v | 0;
                M = M + q + (L >>> 6 ^ L >>> 11 ^ L >>> 25 ^ L << 26 ^ L << 21 ^ L << 7) + (I ^ L & (H ^ I)) + 3329325298 | 0;
                q = I;
                I = H;
                H = L;
                L = x + M | 0;
                x = U;
                U = k;
                k = _;
                _ = M + (k & U ^ x & (k ^ U)) + (k >>> 2 ^ k >>> 13 ^ k >>> 22 ^ k << 30 ^ k << 19 ^ k << 10) | 0;
                i = i + _ | 0;
                n = n + k | 0;
                s = s + U | 0;
                a = a + x | 0;
                h = h + L | 0;
                o = o + H | 0;
                f = f + I | 0;
                u = u + q | 0;
            }
            function I(t) {
                t = t | 0;
                H(L[t | 0] << 24 | L[t | 1] << 16 | L[t | 2] << 8 | L[t | 3], L[t | 4] << 24 | L[t | 5] << 16 | L[t | 6] << 8 | L[t | 7], L[t | 8] << 24 | L[t | 9] << 16 | L[t | 10] << 8 | L[t | 11], L[t | 12] << 24 | L[t | 13] << 16 | L[t | 14] << 8 | L[t | 15], L[t | 16] << 24 | L[t | 17] << 16 | L[t | 18] << 8 | L[t | 19], L[t | 20] << 24 | L[t | 21] << 16 | L[t | 22] << 8 | L[t | 23], L[t | 24] << 24 | L[t | 25] << 16 | L[t | 26] << 8 | L[t | 27], L[t | 28] << 24 | L[t | 29] << 16 | L[t | 30] << 8 | L[t | 31], L[t | 32] << 24 | L[t | 33] << 16 | L[t | 34] << 8 | L[t | 35], L[t | 36] << 24 | L[t | 37] << 16 | L[t | 38] << 8 | L[t | 39], L[t | 40] << 24 | L[t | 41] << 16 | L[t | 42] << 8 | L[t | 43], L[t | 44] << 24 | L[t | 45] << 16 | L[t | 46] << 8 | L[t | 47], L[t | 48] << 24 | L[t | 49] << 16 | L[t | 50] << 8 | L[t | 51], L[t | 52] << 24 | L[t | 53] << 16 | L[t | 54] << 8 | L[t | 55], L[t | 56] << 24 | L[t | 57] << 16 | L[t | 58] << 8 | L[t | 59], L[t | 60] << 24 | L[t | 61] << 16 | L[t | 62] << 8 | L[t | 63]);
            }
            function q(t) {
                t = t | 0;
                L[t | 0] = i >>> 24;
                L[t | 1] = i >>> 16 & 255;
                L[t | 2] = i >>> 8 & 255;
                L[t | 3] = i & 255;
                L[t | 4] = n >>> 24;
                L[t | 5] = n >>> 16 & 255;
                L[t | 6] = n >>> 8 & 255;
                L[t | 7] = n & 255;
                L[t | 8] = s >>> 24;
                L[t | 9] = s >>> 16 & 255;
                L[t | 10] = s >>> 8 & 255;
                L[t | 11] = s & 255;
                L[t | 12] = a >>> 24;
                L[t | 13] = a >>> 16 & 255;
                L[t | 14] = a >>> 8 & 255;
                L[t | 15] = a & 255;
                L[t | 16] = h >>> 24;
                L[t | 17] = h >>> 16 & 255;
                L[t | 18] = h >>> 8 & 255;
                L[t | 19] = h & 255;
                L[t | 20] = o >>> 24;
                L[t | 21] = o >>> 16 & 255;
                L[t | 22] = o >>> 8 & 255;
                L[t | 23] = o & 255;
                L[t | 24] = f >>> 24;
                L[t | 25] = f >>> 16 & 255;
                L[t | 26] = f >>> 8 & 255;
                L[t | 27] = f & 255;
                L[t | 28] = u >>> 24;
                L[t | 29] = u >>> 16 & 255;
                L[t | 30] = u >>> 8 & 255;
                L[t | 31] = u & 255;
            }
            function M() {
                i = 1779033703;
                n = 3144134277;
                s = 1013904242;
                a = 2773480762;
                h = 1359893119;
                o = 2600822924;
                f = 528734635;
                u = 1541459225;
                l = 0;
            }
            function C(t, e, r, c, w, y, p, v, d) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                c = c | 0;
                w = w | 0;
                y = y | 0;
                p = p | 0;
                v = v | 0;
                d = d | 0;
                i = t;
                n = e;
                s = r;
                a = c;
                h = w;
                o = y;
                f = p;
                u = v;
                l = d;
            }
            function Z(t, e) {
                t = t | 0;
                e = e | 0;
                var r = 0;
                if (t & 63) return -1;
                while ((e | 0) >= 64) {
                    I(t);
                    t = t + 64 | 0;
                    e = e - 64 | 0;
                    r = r + 64 | 0;
                }
                l = l + r | 0;
                return r | 0;
            }
            function z(t, e, r) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                var i = 0, n = 0;
                if (t & 63) return -1;
                if (~r) if (r & 31) return -1;
                if ((e | 0) >= 64) {
                    i = Z(t, e) | 0;
                    if ((i | 0) == -1) return -1;
                    t = t + i | 0;
                    e = e - i | 0;
                }
                i = i + e | 0;
                l = l + e | 0;
                L[t | e] = 128;
                if ((e | 0) >= 56) {
                    for (n = e + 1 | 0; (n | 0) < 64; n = n + 1 | 0) L[t | n] = 0;
                    I(t);
                    e = 0;
                    L[t | 0] = 0;
                }
                for (n = e + 1 | 0; (n | 0) < 59; n = n + 1 | 0) L[t | n] = 0;
                L[t | 59] = l >>> 29;
                L[t | 60] = l >>> 21 & 255;
                L[t | 61] = l >>> 13 & 255;
                L[t | 62] = l >>> 5 & 255;
                L[t | 63] = l << 3 & 255;
                I(t);
                if (~r) q(r);
                return i | 0;
            }
            function O() {
                i = c;
                n = w;
                s = y;
                a = p;
                h = v;
                o = d;
                f = g;
                u = b;
                l = 64;
            }
            function T() {
                i = m;
                n = S;
                s = A;
                a = E;
                h = _;
                o = k;
                f = U;
                u = x;
                l = 64;
            }
            function R(t, e, r, L, I, q, C, Z, z, O, T, R, B, P, K, N) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                L = L | 0;
                I = I | 0;
                q = q | 0;
                C = C | 0;
                Z = Z | 0;
                z = z | 0;
                O = O | 0;
                T = T | 0;
                R = R | 0;
                B = B | 0;
                P = P | 0;
                K = K | 0;
                N = N | 0;
                M();
                H(t ^ 1549556828, e ^ 1549556828, r ^ 1549556828, L ^ 1549556828, I ^ 1549556828, q ^ 1549556828, C ^ 1549556828, Z ^ 1549556828, z ^ 1549556828, O ^ 1549556828, T ^ 1549556828, R ^ 1549556828, B ^ 1549556828, P ^ 1549556828, K ^ 1549556828, N ^ 1549556828);
                m = i;
                S = n;
                A = s;
                E = a;
                _ = h;
                k = o;
                U = f;
                x = u;
                M();
                H(t ^ 909522486, e ^ 909522486, r ^ 909522486, L ^ 909522486, I ^ 909522486, q ^ 909522486, C ^ 909522486, Z ^ 909522486, z ^ 909522486, O ^ 909522486, T ^ 909522486, R ^ 909522486, B ^ 909522486, P ^ 909522486, K ^ 909522486, N ^ 909522486);
                c = i;
                w = n;
                y = s;
                p = a;
                v = h;
                d = o;
                g = f;
                b = u;
                l = 64;
            }
            function B(t, e, r) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                var l = 0, c = 0, w = 0, y = 0, p = 0, v = 0, d = 0, g = 0, b = 0;
                if (t & 63) return -1;
                if (~r) if (r & 31) return -1;
                b = z(t, e, -1) | 0;
                l = i, c = n, w = s, y = a, p = h, v = o, d = f, g = u;
                T();
                H(l, c, w, y, p, v, d, g, 2147483648, 0, 0, 0, 0, 0, 0, 768);
                if (~r) q(r);
                return b | 0;
            }
            function P(t, e, r, l, c) {
                t = t | 0;
                e = e | 0;
                r = r | 0;
                l = l | 0;
                c = c | 0;
                var w = 0, y = 0, p = 0, v = 0, d = 0, g = 0, b = 0, m = 0, S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0, I = 0;
                if (t & 63) return -1;
                if (~c) if (c & 31) return -1;
                L[t + e | 0] = r >>> 24;
                L[t + e + 1 | 0] = r >>> 16 & 255;
                L[t + e + 2 | 0] = r >>> 8 & 255;
                L[t + e + 3 | 0] = r & 255;
                B(t, e + 4 | 0, -1) | 0;
                w = S = i, y = A = n, p = E = s, v = _ = a, d = k = h, g = U = o, b = x = f, m = I = u;
                l = l - 1 | 0;
                while ((l | 0) > 0) {
                    O();
                    H(S, A, E, _, k, U, x, I, 2147483648, 0, 0, 0, 0, 0, 0, 768);
                    S = i, A = n, E = s, _ = a, k = h, U = o, x = f, I = u;
                    T();
                    H(S, A, E, _, k, U, x, I, 2147483648, 0, 0, 0, 0, 0, 0, 768);
                    S = i, A = n, E = s, _ = a, k = h, U = o, x = f, I = u;
                    w = w ^ i;
                    y = y ^ n;
                    p = p ^ s;
                    v = v ^ a;
                    d = d ^ h;
                    g = g ^ o;
                    b = b ^ f;
                    m = m ^ u;
                    l = l - 1 | 0;
                }
                i = w;
                n = y;
                s = p;
                a = v;
                h = d;
                o = g;
                f = b;
                u = m;
                if (~c) q(c);
                return 0;
            }
            return {
                reset: M,
                init: C,
                process: Z,
                finish: z,
                hmac_reset: O,
                hmac_init: R,
                hmac_finish: B,
                pbkdf2_generate_block: P
            };
        }(e, null, this.heap.buffer), this.BLOCK_SIZE = Kt, this.HASH_SIZE = Nt, this.reset();
    }
    function T() {
        return null === Dt && (Dt = new O({
            heapSize: 1048576
        })), Dt;
    }
    function R(t) {
        if (void 0 === t) throw new SyntaxError("data required");
        return T().reset().process(t).finish().result;
    }
    function B(t) {
        if (!(t = t || {}).hash) throw new SyntaxError("option 'hash' is required");
        if (!t.hash.HASH_SIZE) throw new SyntaxError("option 'hash' supplied doesn't seem to be a valid hash function");
        return this.hash = t.hash, this.BLOCK_SIZE = this.hash.BLOCK_SIZE, this.HMAC_SIZE = this.hash.HASH_SIZE, 
        this.key = null, this.verify = null, this.result = null, (void 0 !== t.password || void 0 !== t.verify) && this.reset(t), 
        this;
    }
    function P(t, e) {
        if (c(e) && (e = new Uint8Array(e)), l(e) && (e = s(e)), !w(e)) throw new TypeError("password isn't of expected type");
        var r = new Uint8Array(t.BLOCK_SIZE);
        return r.set(e.length > t.BLOCK_SIZE ? t.reset().process(e).finish().result : e), 
        r;
    }
    function K(t) {
        if (c(t) || w(t)) t = new Uint8Array(t); else {
            if (!l(t)) throw new TypeError("verify tag isn't of expected type");
            t = s(t);
        }
        if (t.length !== this.HMAC_SIZE) throw new i("illegal verification tag size");
        this.verify = t;
    }
    function N(t) {
        if (null === this.key) throw new r("no key is associated with the instance");
        if (null !== this.result) throw new r("state must be reset before processing new data");
        return this.hash.process(t), this;
    }
    function j(t) {
        return (t = t || {}).hash instanceof C || (t.hash = Z()), B.call(this, t), this;
    }
    function D() {
        return null === Wt && (Wt = new j()), Wt;
    }
    function G(t) {
        return (t = t || {}).hash instanceof O || (t.hash = T()), B.call(this, t), this;
    }
    function F() {
        return null === Jt && (Jt = new G()), Jt;
    }
    function W(t, e) {
        if (void 0 === t) throw new SyntaxError("data required");
        if (void 0 === e) throw new SyntaxError("password required");
        return D().reset({
            password: e
        }).process(t).finish().result;
    }
    function V(t, e) {
        if (void 0 === t) throw new SyntaxError("data required");
        if (void 0 === e) throw new SyntaxError("password required");
        return F().reset({
            password: e
        }).process(t).finish().result;
    }
    function J(t) {
        if (!(t = t || {}).hmac) throw new SyntaxError("option 'hmac' is required");
        if (!t.hmac.HMAC_SIZE) throw new SyntaxError("option 'hmac' supplied doesn't seem to be a valid HMAC function");
        this.hmac = t.hmac, this.count = t.count || 4096, this.length = t.length || this.hmac.HMAC_SIZE, 
        this.result = null;
        var e = t.password;
        return (e || l(e)) && this.reset(t), this;
    }
    function Q(t) {
        return this.result = null, this.hmac.reset(t), this;
    }
    function X(t) {
        return (t = t || {}).hmac instanceof j || (t.hmac = D()), J.call(this, t), this;
    }
    function Y(t) {
        return (t = t || {}).hmac instanceof G || (t.hmac = F()), J.call(this, t), this;
    }
    function $() {
        return null === te && (te = new Y()), te;
    }
    function tt(t, e, r, i) {
        if (void 0 === t) throw new SyntaxError("password required");
        if (void 0 === e) throw new SyntaxError("salt required");
        return (null === Yt && (Yt = new X()), Yt).reset({
            password: t
        }).generate(e, r, i).result;
    }
    function et(t, e, r, i) {
        if (void 0 === t) throw new SyntaxError("password required");
        if (void 0 === e) throw new SyntaxError("salt required");
        return $().reset({
            password: t
        }).generate(e, r, i).result;
    }
    function rt() {
        if (void 0 !== he) i = new Uint8Array(32), ee.call(he, i), le(i); else {
            var t, r, i = new Et(3);
            i[0] = se(), i[1] = ne(), i[2] = oe(), i = new Uint8Array(i.buffer);
            var n = $();
            for (t = 0; 100 > t; t++) i = n.reset({
                password: i
            }).generate(e.location.href, 1e3, 32).result, r = oe(), i[0] ^= r >>> 24, i[1] ^= r >>> 16, 
            i[2] ^= r >>> 8, i[3] ^= r;
            le(i);
        }
        ce = 0, we = !0;
    }
    function it(t) {
        if (!c(t) && !y(t)) throw new TypeError("bad seed type");
        var e = t.byteOffest || 0, r = t.byteLength || t.length, i = new Uint8Array(t.buffer || t, e, r);
        le(i), ce = 0;
        for (var n = 0, s = 0; s < i.length; s++) n |= i[s], i[s] = 0;
        return 0 !== n && (pe += 4 * r), ye = pe >= ve;
    }
    function nt(t) {
        if (we || rt(), !ye && void 0 === he) {
            if (!de) throw new n("No strong PRNGs available. Use asmCrypto.random.seed().");
            void 0 !== ie && ie.error("No strong PRNGs available; your security is greatly lowered. Use asmCrypto.random.seed().");
        }
        if (!ge && !ye && void 0 !== he && void 0 !== ie) {
            var e = new Error().stack;
            be[e] |= 0, be[e]++ || ie.warn("asmCrypto PRNG not seeded; your security relies on your system PRNG. If this is not acceptable, use asmCrypto.random.seed().");
        }
        if (!c(t) && !y(t)) throw new TypeError("unexpected buffer type");
        var r, i, s = t.byteOffset || 0, a = t.byteLength || t.length, h = new Uint8Array(t.buffer || t, s, a);
        for (void 0 !== he && ee.call(he, h), r = 0; a > r; r++) 0 == (3 & r) && (ce >= 1099511627776 && rt(), 
        i = ue(), ce++), h[r] ^= i, i >>>= 8;
    }
    function st() {
        (!we || ce >= 1099511627776) && rt();
        var t = (1048576 * ue() + (ue() >>> 12)) / 4503599627370496;
        return ce += 2, t;
    }
    function at(t, e) {
        return t * e | 0;
    }
    function ht(t, e, r) {
        "use asm";
        var i = 0;
        var n = new t.Uint32Array(r);
        var s = t.Math.imul;
        function a(t) {
            t = t | 0;
            i = t = t + 31 & -32;
            return t | 0;
        }
        function h(t) {
            t = t | 0;
            var e = 0;
            e = i;
            i = e + (t + 31 & -32) | 0;
            return e | 0;
        }
        function o(t) {
            t = t | 0;
            i = i - (t + 31 & -32) | 0;
        }
        function f(t, e, r) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            var i = 0;
            if ((e | 0) > (r | 0)) {
                for (;(i | 0) < (t | 0); i = i + 4 | 0) {
                    n[r + i >> 2] = n[e + i >> 2];
                }
            } else {
                for (i = t - 4 | 0; (i | 0) >= 0; i = i - 4 | 0) {
                    n[r + i >> 2] = n[e + i >> 2];
                }
            }
        }
        function u(t, e, r) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            var i = 0;
            for (;(i | 0) < (t | 0); i = i + 4 | 0) {
                n[r + i >> 2] = e;
            }
        }
        function l(t, e, r, i) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            var s = 0, a = 0, h = 0, o = 0, f = 0;
            if ((i | 0) <= 0) i = e;
            if ((i | 0) < (e | 0)) e = i;
            a = 1;
            for (;(f | 0) < (e | 0); f = f + 4 | 0) {
                s = ~n[t + f >> 2];
                h = (s & 65535) + a | 0;
                o = (s >>> 16) + (h >>> 16) | 0;
                n[r + f >> 2] = o << 16 | h & 65535;
                a = o >>> 16;
            }
            for (;(f | 0) < (i | 0); f = f + 4 | 0) {
                n[r + f >> 2] = a - 1 | 0;
            }
            return a | 0;
        }
        function c(t, e, r, i) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            var s = 0, a = 0, h = 0;
            if ((e | 0) > (i | 0)) {
                for (h = e - 4 | 0; (h | 0) >= (i | 0); h = h - 4 | 0) {
                    if (n[t + h >> 2] | 0) return 1;
                }
            } else {
                for (h = i - 4 | 0; (h | 0) >= (e | 0); h = h - 4 | 0) {
                    if (n[r + h >> 2] | 0) return -1;
                }
            }
            for (;(h | 0) >= 0; h = h - 4 | 0) {
                s = n[t + h >> 2] | 0, a = n[r + h >> 2] | 0;
                if (s >>> 0 < a >>> 0) return -1;
                if (s >>> 0 > a >>> 0) return 1;
            }
            return 0;
        }
        function w(t, e) {
            t = t | 0;
            e = e | 0;
            var r = 0;
            for (r = e - 4 | 0; (r | 0) >= 0; r = r - 4 | 0) {
                if (n[t + r >> 2] | 0) return r + 4 | 0;
            }
            return 0;
        }
        function y(t, e, r, i, s, a) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            s = s | 0;
            a = a | 0;
            var h = 0, o = 0, f = 0, u = 0, l = 0, c = 0;
            if ((e | 0) < (i | 0)) {
                u = t, t = r, r = u;
                u = e, e = i, i = u;
            }
            if ((a | 0) <= 0) a = e + 4 | 0;
            if ((a | 0) < (i | 0)) e = i = a;
            for (;(c | 0) < (i | 0); c = c + 4 | 0) {
                h = n[t + c >> 2] | 0;
                o = n[r + c >> 2] | 0;
                u = ((h & 65535) + (o & 65535) | 0) + f | 0;
                l = ((h >>> 16) + (o >>> 16) | 0) + (u >>> 16) | 0;
                n[s + c >> 2] = u & 65535 | l << 16;
                f = l >>> 16;
            }
            for (;(c | 0) < (e | 0); c = c + 4 | 0) {
                h = n[t + c >> 2] | 0;
                u = (h & 65535) + f | 0;
                l = (h >>> 16) + (u >>> 16) | 0;
                n[s + c >> 2] = u & 65535 | l << 16;
                f = l >>> 16;
            }
            for (;(c | 0) < (a | 0); c = c + 4 | 0) {
                n[s + c >> 2] = f | 0;
                f = 0;
            }
            return f | 0;
        }
        function p(t, e, r, i, s, a) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            s = s | 0;
            a = a | 0;
            var h = 0, o = 0, f = 0, u = 0, l = 0, c = 0;
            if ((a | 0) <= 0) a = (e | 0) > (i | 0) ? e + 4 | 0 : i + 4 | 0;
            if ((a | 0) < (e | 0)) e = a;
            if ((a | 0) < (i | 0)) i = a;
            if ((e | 0) < (i | 0)) {
                for (;(c | 0) < (e | 0); c = c + 4 | 0) {
                    h = n[t + c >> 2] | 0;
                    o = n[r + c >> 2] | 0;
                    u = ((h & 65535) - (o & 65535) | 0) + f | 0;
                    l = ((h >>> 16) - (o >>> 16) | 0) + (u >> 16) | 0;
                    n[s + c >> 2] = u & 65535 | l << 16;
                    f = l >> 16;
                }
                for (;(c | 0) < (i | 0); c = c + 4 | 0) {
                    o = n[r + c >> 2] | 0;
                    u = f - (o & 65535) | 0;
                    l = (u >> 16) - (o >>> 16) | 0;
                    n[s + c >> 2] = u & 65535 | l << 16;
                    f = l >> 16;
                }
            } else {
                for (;(c | 0) < (i | 0); c = c + 4 | 0) {
                    h = n[t + c >> 2] | 0;
                    o = n[r + c >> 2] | 0;
                    u = ((h & 65535) - (o & 65535) | 0) + f | 0;
                    l = ((h >>> 16) - (o >>> 16) | 0) + (u >> 16) | 0;
                    n[s + c >> 2] = u & 65535 | l << 16;
                    f = l >> 16;
                }
                for (;(c | 0) < (e | 0); c = c + 4 | 0) {
                    h = n[t + c >> 2] | 0;
                    u = (h & 65535) + f | 0;
                    l = (h >>> 16) + (u >> 16) | 0;
                    n[s + c >> 2] = u & 65535 | l << 16;
                    f = l >> 16;
                }
            }
            for (;(c | 0) < (a | 0); c = c + 4 | 0) {
                n[s + c >> 2] = f | 0;
            }
            return f | 0;
        }
        function v(t, e, r, i, a, h) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            a = a | 0;
            h = h | 0;
            var o = 0, f = 0, u = 0, l = 0, c = 0, w = 0, y = 0, p = 0, v = 0, d = 0, g = 0, b = 0, m = 0, S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0, L = 0, H = 0, I = 0, q = 0, M = 0, C = 0, Z = 0, z = 0, O = 0, T = 0, R = 0, B = 0, P = 0, K = 0, N = 0, j = 0, D = 0, G = 0, F = 0, W = 0, V = 0, J = 0, Q = 0, X = 0, Y = 0, $ = 0, tt = 0, et = 0, rt = 0, it = 0, nt = 0, st = 0, at = 0, ht = 0, ot = 0, ft = 0, ut = 0;
            if ((e | 0) > (i | 0)) {
                rt = t, it = e;
                t = r, e = i;
                r = rt, i = it;
            }
            st = e + i | 0;
            if ((h | 0) > (st | 0) | (h | 0) <= 0) h = st;
            if ((h | 0) < (e | 0)) e = h;
            if ((h | 0) < (i | 0)) i = h;
            for (;(at | 0) < (e | 0); at = at + 32 | 0) {
                ht = t + at | 0;
                v = n[(ht | 0) >> 2] | 0, d = n[(ht | 4) >> 2] | 0, g = n[(ht | 8) >> 2] | 0, b = n[(ht | 12) >> 2] | 0, 
                m = n[(ht | 16) >> 2] | 0, S = n[(ht | 20) >> 2] | 0, A = n[(ht | 24) >> 2] | 0, 
                E = n[(ht | 28) >> 2] | 0, o = v & 65535, f = d & 65535, u = g & 65535, l = b & 65535, 
                c = m & 65535, w = S & 65535, y = A & 65535, p = E & 65535, v = v >>> 16, d = d >>> 16, 
                g = g >>> 16, b = b >>> 16, m = m >>> 16, S = S >>> 16, A = A >>> 16, E = E >>> 16;
                V = J = Q = X = Y = $ = tt = et = 0;
                for (ot = 0; (ot | 0) < (i | 0); ot = ot + 32 | 0) {
                    ft = r + ot | 0;
                    ut = a + (at + ot | 0) | 0;
                    M = n[(ft | 0) >> 2] | 0, C = n[(ft | 4) >> 2] | 0, Z = n[(ft | 8) >> 2] | 0, z = n[(ft | 12) >> 2] | 0, 
                    O = n[(ft | 16) >> 2] | 0, T = n[(ft | 20) >> 2] | 0, R = n[(ft | 24) >> 2] | 0, 
                    B = n[(ft | 28) >> 2] | 0, _ = M & 65535, k = C & 65535, U = Z & 65535, x = z & 65535, 
                    L = O & 65535, H = T & 65535, I = R & 65535, q = B & 65535, M = M >>> 16, C = C >>> 16, 
                    Z = Z >>> 16, z = z >>> 16, O = O >>> 16, T = T >>> 16, R = R >>> 16, B = B >>> 16;
                    P = n[(ut | 0) >> 2] | 0, K = n[(ut | 4) >> 2] | 0, N = n[(ut | 8) >> 2] | 0, j = n[(ut | 12) >> 2] | 0, 
                    D = n[(ut | 16) >> 2] | 0, G = n[(ut | 20) >> 2] | 0, F = n[(ut | 24) >> 2] | 0, 
                    W = n[(ut | 28) >> 2] | 0;
                    rt = ((s(o, _) | 0) + (V & 65535) | 0) + (P & 65535) | 0;
                    it = ((s(v, _) | 0) + (V >>> 16) | 0) + (P >>> 16) | 0;
                    nt = ((s(o, M) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(v, M) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    P = nt << 16 | rt & 65535;
                    rt = ((s(o, k) | 0) + (st & 65535) | 0) + (K & 65535) | 0;
                    it = ((s(v, k) | 0) + (st >>> 16) | 0) + (K >>> 16) | 0;
                    nt = ((s(o, C) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(v, C) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    K = nt << 16 | rt & 65535;
                    rt = ((s(o, U) | 0) + (st & 65535) | 0) + (N & 65535) | 0;
                    it = ((s(v, U) | 0) + (st >>> 16) | 0) + (N >>> 16) | 0;
                    nt = ((s(o, Z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(v, Z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    N = nt << 16 | rt & 65535;
                    rt = ((s(o, x) | 0) + (st & 65535) | 0) + (j & 65535) | 0;
                    it = ((s(v, x) | 0) + (st >>> 16) | 0) + (j >>> 16) | 0;
                    nt = ((s(o, z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(v, z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    j = nt << 16 | rt & 65535;
                    rt = ((s(o, L) | 0) + (st & 65535) | 0) + (D & 65535) | 0;
                    it = ((s(v, L) | 0) + (st >>> 16) | 0) + (D >>> 16) | 0;
                    nt = ((s(o, O) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(v, O) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    D = nt << 16 | rt & 65535;
                    rt = ((s(o, H) | 0) + (st & 65535) | 0) + (G & 65535) | 0;
                    it = ((s(v, H) | 0) + (st >>> 16) | 0) + (G >>> 16) | 0;
                    nt = ((s(o, T) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(v, T) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    G = nt << 16 | rt & 65535;
                    rt = ((s(o, I) | 0) + (st & 65535) | 0) + (F & 65535) | 0;
                    it = ((s(v, I) | 0) + (st >>> 16) | 0) + (F >>> 16) | 0;
                    nt = ((s(o, R) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(v, R) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    F = nt << 16 | rt & 65535;
                    rt = ((s(o, q) | 0) + (st & 65535) | 0) + (W & 65535) | 0;
                    it = ((s(v, q) | 0) + (st >>> 16) | 0) + (W >>> 16) | 0;
                    nt = ((s(o, B) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(v, B) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    W = nt << 16 | rt & 65535;
                    V = st;
                    rt = ((s(f, _) | 0) + (J & 65535) | 0) + (K & 65535) | 0;
                    it = ((s(d, _) | 0) + (J >>> 16) | 0) + (K >>> 16) | 0;
                    nt = ((s(f, M) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(d, M) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    K = nt << 16 | rt & 65535;
                    rt = ((s(f, k) | 0) + (st & 65535) | 0) + (N & 65535) | 0;
                    it = ((s(d, k) | 0) + (st >>> 16) | 0) + (N >>> 16) | 0;
                    nt = ((s(f, C) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(d, C) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    N = nt << 16 | rt & 65535;
                    rt = ((s(f, U) | 0) + (st & 65535) | 0) + (j & 65535) | 0;
                    it = ((s(d, U) | 0) + (st >>> 16) | 0) + (j >>> 16) | 0;
                    nt = ((s(f, Z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(d, Z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    j = nt << 16 | rt & 65535;
                    rt = ((s(f, x) | 0) + (st & 65535) | 0) + (D & 65535) | 0;
                    it = ((s(d, x) | 0) + (st >>> 16) | 0) + (D >>> 16) | 0;
                    nt = ((s(f, z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(d, z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    D = nt << 16 | rt & 65535;
                    rt = ((s(f, L) | 0) + (st & 65535) | 0) + (G & 65535) | 0;
                    it = ((s(d, L) | 0) + (st >>> 16) | 0) + (G >>> 16) | 0;
                    nt = ((s(f, O) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(d, O) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    G = nt << 16 | rt & 65535;
                    rt = ((s(f, H) | 0) + (st & 65535) | 0) + (F & 65535) | 0;
                    it = ((s(d, H) | 0) + (st >>> 16) | 0) + (F >>> 16) | 0;
                    nt = ((s(f, T) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(d, T) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    F = nt << 16 | rt & 65535;
                    rt = ((s(f, I) | 0) + (st & 65535) | 0) + (W & 65535) | 0;
                    it = ((s(d, I) | 0) + (st >>> 16) | 0) + (W >>> 16) | 0;
                    nt = ((s(f, R) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(d, R) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    W = nt << 16 | rt & 65535;
                    rt = ((s(f, q) | 0) + (st & 65535) | 0) + (V & 65535) | 0;
                    it = ((s(d, q) | 0) + (st >>> 16) | 0) + (V >>> 16) | 0;
                    nt = ((s(f, B) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(d, B) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    V = nt << 16 | rt & 65535;
                    J = st;
                    rt = ((s(u, _) | 0) + (Q & 65535) | 0) + (N & 65535) | 0;
                    it = ((s(g, _) | 0) + (Q >>> 16) | 0) + (N >>> 16) | 0;
                    nt = ((s(u, M) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(g, M) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    N = nt << 16 | rt & 65535;
                    rt = ((s(u, k) | 0) + (st & 65535) | 0) + (j & 65535) | 0;
                    it = ((s(g, k) | 0) + (st >>> 16) | 0) + (j >>> 16) | 0;
                    nt = ((s(u, C) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(g, C) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    j = nt << 16 | rt & 65535;
                    rt = ((s(u, U) | 0) + (st & 65535) | 0) + (D & 65535) | 0;
                    it = ((s(g, U) | 0) + (st >>> 16) | 0) + (D >>> 16) | 0;
                    nt = ((s(u, Z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(g, Z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    D = nt << 16 | rt & 65535;
                    rt = ((s(u, x) | 0) + (st & 65535) | 0) + (G & 65535) | 0;
                    it = ((s(g, x) | 0) + (st >>> 16) | 0) + (G >>> 16) | 0;
                    nt = ((s(u, z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(g, z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    G = nt << 16 | rt & 65535;
                    rt = ((s(u, L) | 0) + (st & 65535) | 0) + (F & 65535) | 0;
                    it = ((s(g, L) | 0) + (st >>> 16) | 0) + (F >>> 16) | 0;
                    nt = ((s(u, O) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(g, O) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    F = nt << 16 | rt & 65535;
                    rt = ((s(u, H) | 0) + (st & 65535) | 0) + (W & 65535) | 0;
                    it = ((s(g, H) | 0) + (st >>> 16) | 0) + (W >>> 16) | 0;
                    nt = ((s(u, T) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(g, T) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    W = nt << 16 | rt & 65535;
                    rt = ((s(u, I) | 0) + (st & 65535) | 0) + (V & 65535) | 0;
                    it = ((s(g, I) | 0) + (st >>> 16) | 0) + (V >>> 16) | 0;
                    nt = ((s(u, R) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(g, R) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    V = nt << 16 | rt & 65535;
                    rt = ((s(u, q) | 0) + (st & 65535) | 0) + (J & 65535) | 0;
                    it = ((s(g, q) | 0) + (st >>> 16) | 0) + (J >>> 16) | 0;
                    nt = ((s(u, B) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(g, B) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    J = nt << 16 | rt & 65535;
                    Q = st;
                    rt = ((s(l, _) | 0) + (X & 65535) | 0) + (j & 65535) | 0;
                    it = ((s(b, _) | 0) + (X >>> 16) | 0) + (j >>> 16) | 0;
                    nt = ((s(l, M) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(b, M) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    j = nt << 16 | rt & 65535;
                    rt = ((s(l, k) | 0) + (st & 65535) | 0) + (D & 65535) | 0;
                    it = ((s(b, k) | 0) + (st >>> 16) | 0) + (D >>> 16) | 0;
                    nt = ((s(l, C) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(b, C) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    D = nt << 16 | rt & 65535;
                    rt = ((s(l, U) | 0) + (st & 65535) | 0) + (G & 65535) | 0;
                    it = ((s(b, U) | 0) + (st >>> 16) | 0) + (G >>> 16) | 0;
                    nt = ((s(l, Z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(b, Z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    G = nt << 16 | rt & 65535;
                    rt = ((s(l, x) | 0) + (st & 65535) | 0) + (F & 65535) | 0;
                    it = ((s(b, x) | 0) + (st >>> 16) | 0) + (F >>> 16) | 0;
                    nt = ((s(l, z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(b, z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    F = nt << 16 | rt & 65535;
                    rt = ((s(l, L) | 0) + (st & 65535) | 0) + (W & 65535) | 0;
                    it = ((s(b, L) | 0) + (st >>> 16) | 0) + (W >>> 16) | 0;
                    nt = ((s(l, O) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(b, O) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    W = nt << 16 | rt & 65535;
                    rt = ((s(l, H) | 0) + (st & 65535) | 0) + (V & 65535) | 0;
                    it = ((s(b, H) | 0) + (st >>> 16) | 0) + (V >>> 16) | 0;
                    nt = ((s(l, T) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(b, T) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    V = nt << 16 | rt & 65535;
                    rt = ((s(l, I) | 0) + (st & 65535) | 0) + (J & 65535) | 0;
                    it = ((s(b, I) | 0) + (st >>> 16) | 0) + (J >>> 16) | 0;
                    nt = ((s(l, R) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(b, R) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    J = nt << 16 | rt & 65535;
                    rt = ((s(l, q) | 0) + (st & 65535) | 0) + (Q & 65535) | 0;
                    it = ((s(b, q) | 0) + (st >>> 16) | 0) + (Q >>> 16) | 0;
                    nt = ((s(l, B) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(b, B) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    Q = nt << 16 | rt & 65535;
                    X = st;
                    rt = ((s(c, _) | 0) + (Y & 65535) | 0) + (D & 65535) | 0;
                    it = ((s(m, _) | 0) + (Y >>> 16) | 0) + (D >>> 16) | 0;
                    nt = ((s(c, M) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(m, M) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    D = nt << 16 | rt & 65535;
                    rt = ((s(c, k) | 0) + (st & 65535) | 0) + (G & 65535) | 0;
                    it = ((s(m, k) | 0) + (st >>> 16) | 0) + (G >>> 16) | 0;
                    nt = ((s(c, C) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(m, C) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    G = nt << 16 | rt & 65535;
                    rt = ((s(c, U) | 0) + (st & 65535) | 0) + (F & 65535) | 0;
                    it = ((s(m, U) | 0) + (st >>> 16) | 0) + (F >>> 16) | 0;
                    nt = ((s(c, Z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(m, Z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    F = nt << 16 | rt & 65535;
                    rt = ((s(c, x) | 0) + (st & 65535) | 0) + (W & 65535) | 0;
                    it = ((s(m, x) | 0) + (st >>> 16) | 0) + (W >>> 16) | 0;
                    nt = ((s(c, z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(m, z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    W = nt << 16 | rt & 65535;
                    rt = ((s(c, L) | 0) + (st & 65535) | 0) + (V & 65535) | 0;
                    it = ((s(m, L) | 0) + (st >>> 16) | 0) + (V >>> 16) | 0;
                    nt = ((s(c, O) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(m, O) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    V = nt << 16 | rt & 65535;
                    rt = ((s(c, H) | 0) + (st & 65535) | 0) + (J & 65535) | 0;
                    it = ((s(m, H) | 0) + (st >>> 16) | 0) + (J >>> 16) | 0;
                    nt = ((s(c, T) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(m, T) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    J = nt << 16 | rt & 65535;
                    rt = ((s(c, I) | 0) + (st & 65535) | 0) + (Q & 65535) | 0;
                    it = ((s(m, I) | 0) + (st >>> 16) | 0) + (Q >>> 16) | 0;
                    nt = ((s(c, R) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(m, R) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    Q = nt << 16 | rt & 65535;
                    rt = ((s(c, q) | 0) + (st & 65535) | 0) + (X & 65535) | 0;
                    it = ((s(m, q) | 0) + (st >>> 16) | 0) + (X >>> 16) | 0;
                    nt = ((s(c, B) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(m, B) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    X = nt << 16 | rt & 65535;
                    Y = st;
                    rt = ((s(w, _) | 0) + ($ & 65535) | 0) + (G & 65535) | 0;
                    it = ((s(S, _) | 0) + ($ >>> 16) | 0) + (G >>> 16) | 0;
                    nt = ((s(w, M) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(S, M) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    G = nt << 16 | rt & 65535;
                    rt = ((s(w, k) | 0) + (st & 65535) | 0) + (F & 65535) | 0;
                    it = ((s(S, k) | 0) + (st >>> 16) | 0) + (F >>> 16) | 0;
                    nt = ((s(w, C) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(S, C) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    F = nt << 16 | rt & 65535;
                    rt = ((s(w, U) | 0) + (st & 65535) | 0) + (W & 65535) | 0;
                    it = ((s(S, U) | 0) + (st >>> 16) | 0) + (W >>> 16) | 0;
                    nt = ((s(w, Z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(S, Z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    W = nt << 16 | rt & 65535;
                    rt = ((s(w, x) | 0) + (st & 65535) | 0) + (V & 65535) | 0;
                    it = ((s(S, x) | 0) + (st >>> 16) | 0) + (V >>> 16) | 0;
                    nt = ((s(w, z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(S, z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    V = nt << 16 | rt & 65535;
                    rt = ((s(w, L) | 0) + (st & 65535) | 0) + (J & 65535) | 0;
                    it = ((s(S, L) | 0) + (st >>> 16) | 0) + (J >>> 16) | 0;
                    nt = ((s(w, O) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(S, O) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    J = nt << 16 | rt & 65535;
                    rt = ((s(w, H) | 0) + (st & 65535) | 0) + (Q & 65535) | 0;
                    it = ((s(S, H) | 0) + (st >>> 16) | 0) + (Q >>> 16) | 0;
                    nt = ((s(w, T) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(S, T) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    Q = nt << 16 | rt & 65535;
                    rt = ((s(w, I) | 0) + (st & 65535) | 0) + (X & 65535) | 0;
                    it = ((s(S, I) | 0) + (st >>> 16) | 0) + (X >>> 16) | 0;
                    nt = ((s(w, R) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(S, R) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    X = nt << 16 | rt & 65535;
                    rt = ((s(w, q) | 0) + (st & 65535) | 0) + (Y & 65535) | 0;
                    it = ((s(S, q) | 0) + (st >>> 16) | 0) + (Y >>> 16) | 0;
                    nt = ((s(w, B) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(S, B) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    Y = nt << 16 | rt & 65535;
                    $ = st;
                    rt = ((s(y, _) | 0) + (tt & 65535) | 0) + (F & 65535) | 0;
                    it = ((s(A, _) | 0) + (tt >>> 16) | 0) + (F >>> 16) | 0;
                    nt = ((s(y, M) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(A, M) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    F = nt << 16 | rt & 65535;
                    rt = ((s(y, k) | 0) + (st & 65535) | 0) + (W & 65535) | 0;
                    it = ((s(A, k) | 0) + (st >>> 16) | 0) + (W >>> 16) | 0;
                    nt = ((s(y, C) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(A, C) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    W = nt << 16 | rt & 65535;
                    rt = ((s(y, U) | 0) + (st & 65535) | 0) + (V & 65535) | 0;
                    it = ((s(A, U) | 0) + (st >>> 16) | 0) + (V >>> 16) | 0;
                    nt = ((s(y, Z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(A, Z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    V = nt << 16 | rt & 65535;
                    rt = ((s(y, x) | 0) + (st & 65535) | 0) + (J & 65535) | 0;
                    it = ((s(A, x) | 0) + (st >>> 16) | 0) + (J >>> 16) | 0;
                    nt = ((s(y, z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(A, z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    J = nt << 16 | rt & 65535;
                    rt = ((s(y, L) | 0) + (st & 65535) | 0) + (Q & 65535) | 0;
                    it = ((s(A, L) | 0) + (st >>> 16) | 0) + (Q >>> 16) | 0;
                    nt = ((s(y, O) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(A, O) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    Q = nt << 16 | rt & 65535;
                    rt = ((s(y, H) | 0) + (st & 65535) | 0) + (X & 65535) | 0;
                    it = ((s(A, H) | 0) + (st >>> 16) | 0) + (X >>> 16) | 0;
                    nt = ((s(y, T) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(A, T) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    X = nt << 16 | rt & 65535;
                    rt = ((s(y, I) | 0) + (st & 65535) | 0) + (Y & 65535) | 0;
                    it = ((s(A, I) | 0) + (st >>> 16) | 0) + (Y >>> 16) | 0;
                    nt = ((s(y, R) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(A, R) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    Y = nt << 16 | rt & 65535;
                    rt = ((s(y, q) | 0) + (st & 65535) | 0) + ($ & 65535) | 0;
                    it = ((s(A, q) | 0) + (st >>> 16) | 0) + ($ >>> 16) | 0;
                    nt = ((s(y, B) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(A, B) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    $ = nt << 16 | rt & 65535;
                    tt = st;
                    rt = ((s(p, _) | 0) + (et & 65535) | 0) + (W & 65535) | 0;
                    it = ((s(E, _) | 0) + (et >>> 16) | 0) + (W >>> 16) | 0;
                    nt = ((s(p, M) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(E, M) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    W = nt << 16 | rt & 65535;
                    rt = ((s(p, k) | 0) + (st & 65535) | 0) + (V & 65535) | 0;
                    it = ((s(E, k) | 0) + (st >>> 16) | 0) + (V >>> 16) | 0;
                    nt = ((s(p, C) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(E, C) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    V = nt << 16 | rt & 65535;
                    rt = ((s(p, U) | 0) + (st & 65535) | 0) + (J & 65535) | 0;
                    it = ((s(E, U) | 0) + (st >>> 16) | 0) + (J >>> 16) | 0;
                    nt = ((s(p, Z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(E, Z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    J = nt << 16 | rt & 65535;
                    rt = ((s(p, x) | 0) + (st & 65535) | 0) + (Q & 65535) | 0;
                    it = ((s(E, x) | 0) + (st >>> 16) | 0) + (Q >>> 16) | 0;
                    nt = ((s(p, z) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(E, z) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    Q = nt << 16 | rt & 65535;
                    rt = ((s(p, L) | 0) + (st & 65535) | 0) + (X & 65535) | 0;
                    it = ((s(E, L) | 0) + (st >>> 16) | 0) + (X >>> 16) | 0;
                    nt = ((s(p, O) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(E, O) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    X = nt << 16 | rt & 65535;
                    rt = ((s(p, H) | 0) + (st & 65535) | 0) + (Y & 65535) | 0;
                    it = ((s(E, H) | 0) + (st >>> 16) | 0) + (Y >>> 16) | 0;
                    nt = ((s(p, T) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(E, T) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    Y = nt << 16 | rt & 65535;
                    rt = ((s(p, I) | 0) + (st & 65535) | 0) + ($ & 65535) | 0;
                    it = ((s(E, I) | 0) + (st >>> 16) | 0) + ($ >>> 16) | 0;
                    nt = ((s(p, R) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(E, R) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    $ = nt << 16 | rt & 65535;
                    rt = ((s(p, q) | 0) + (st & 65535) | 0) + (tt & 65535) | 0;
                    it = ((s(E, q) | 0) + (st >>> 16) | 0) + (tt >>> 16) | 0;
                    nt = ((s(p, B) | 0) + (it & 65535) | 0) + (rt >>> 16) | 0;
                    st = ((s(E, B) | 0) + (it >>> 16) | 0) + (nt >>> 16) | 0;
                    tt = nt << 16 | rt & 65535;
                    et = st;
                    n[(ut | 0) >> 2] = P, n[(ut | 4) >> 2] = K, n[(ut | 8) >> 2] = N, n[(ut | 12) >> 2] = j, 
                    n[(ut | 16) >> 2] = D, n[(ut | 20) >> 2] = G, n[(ut | 24) >> 2] = F, n[(ut | 28) >> 2] = W;
                }
                ut = a + (at + ot | 0) | 0;
                n[(ut | 0) >> 2] = V, n[(ut | 4) >> 2] = J, n[(ut | 8) >> 2] = Q, n[(ut | 12) >> 2] = X, 
                n[(ut | 16) >> 2] = Y, n[(ut | 20) >> 2] = $, n[(ut | 24) >> 2] = tt, n[(ut | 28) >> 2] = et;
            }
        }
        function d(t, e, r) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            var i = 0, a = 0, h = 0, o = 0, f = 0, u = 0, l = 0, c = 0, w = 0, y = 0, p = 0, v = 0, d = 0, g = 0, b = 0, m = 0, S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0, L = 0, H = 0, I = 0, q = 0, M = 0, C = 0, Z = 0, z = 0, O = 0, T = 0, R = 0, B = 0, P = 0, K = 0, N = 0, j = 0, D = 0, G = 0, F = 0, W = 0, V = 0, J = 0, Q = 0, X = 0, Y = 0, $ = 0, tt = 0, et = 0, rt = 0, it = 0, nt = 0, st = 0, at = 0, ht = 0, ot = 0, ft = 0, ut = 0, lt = 0, ct = 0, wt = 0, yt = 0;
            for (;(ft | 0) < (e | 0); ft = ft + 4 | 0) {
                yt = r + (ft << 1) | 0;
                w = n[t + ft >> 2] | 0, i = w & 65535, w = w >>> 16;
                $ = s(i, i) | 0;
                tt = (s(i, w) | 0) + ($ >>> 17) | 0;
                et = (s(w, w) | 0) + (tt >>> 15) | 0;
                n[yt >> 2] = tt << 17 | $ & 131071;
                n[(yt | 4) >> 2] = et;
            }
            for (ot = 0; (ot | 0) < (e | 0); ot = ot + 8 | 0) {
                ct = t + ot | 0, yt = r + (ot << 1) | 0;
                w = n[ct >> 2] | 0, i = w & 65535, w = w >>> 16;
                H = n[(ct | 4) >> 2] | 0, S = H & 65535, H = H >>> 16;
                $ = s(i, S) | 0;
                tt = (s(i, H) | 0) + ($ >>> 16) | 0;
                et = (s(w, S) | 0) + (tt & 65535) | 0;
                nt = ((s(w, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                st = n[(yt | 4) >> 2] | 0;
                $ = (st & 65535) + (($ & 65535) << 1) | 0;
                et = ((st >>> 16) + ((et & 65535) << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 4) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[(yt | 8) >> 2] | 0;
                $ = ((st & 65535) + ((nt & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (nt >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 8) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                if (rt) {
                    st = n[(yt | 12) >> 2] | 0;
                    $ = (st & 65535) + rt | 0;
                    et = (st >>> 16) + ($ >>> 16) | 0;
                    n[(yt | 12) >> 2] = et << 16 | $ & 65535;
                }
            }
            for (ot = 0; (ot | 0) < (e | 0); ot = ot + 16 | 0) {
                ct = t + ot | 0, yt = r + (ot << 1) | 0;
                w = n[ct >> 2] | 0, i = w & 65535, w = w >>> 16, y = n[(ct | 4) >> 2] | 0, a = y & 65535, 
                y = y >>> 16;
                H = n[(ct | 8) >> 2] | 0, S = H & 65535, H = H >>> 16, I = n[(ct | 12) >> 2] | 0, 
                A = I & 65535, I = I >>> 16;
                $ = s(i, S) | 0;
                tt = s(w, S) | 0;
                et = ((s(i, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(w, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                T = et << 16 | $ & 65535;
                $ = (s(i, A) | 0) + (nt & 65535) | 0;
                tt = (s(w, A) | 0) + (nt >>> 16) | 0;
                et = ((s(i, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(w, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                R = et << 16 | $ & 65535;
                B = nt;
                $ = (s(a, S) | 0) + (R & 65535) | 0;
                tt = (s(y, S) | 0) + (R >>> 16) | 0;
                et = ((s(a, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(y, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                R = et << 16 | $ & 65535;
                $ = ((s(a, A) | 0) + (B & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(y, A) | 0) + (B >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(a, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(y, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                B = et << 16 | $ & 65535;
                P = nt;
                st = n[(yt | 8) >> 2] | 0;
                $ = (st & 65535) + ((T & 65535) << 1) | 0;
                et = ((st >>> 16) + (T >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 8) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[(yt | 12) >> 2] | 0;
                $ = ((st & 65535) + ((R & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (R >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 12) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[(yt | 16) >> 2] | 0;
                $ = ((st & 65535) + ((B & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (B >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 16) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[(yt | 20) >> 2] | 0;
                $ = ((st & 65535) + ((P & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (P >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 20) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                for (lt = 24; !!rt & (lt | 0) < 32; lt = lt + 4 | 0) {
                    st = n[(yt | lt) >> 2] | 0;
                    $ = (st & 65535) + rt | 0;
                    et = (st >>> 16) + ($ >>> 16) | 0;
                    n[(yt | lt) >> 2] = et << 16 | $ & 65535;
                    rt = et >>> 16;
                }
            }
            for (ot = 0; (ot | 0) < (e | 0); ot = ot + 32 | 0) {
                ct = t + ot | 0, yt = r + (ot << 1) | 0;
                w = n[ct >> 2] | 0, i = w & 65535, w = w >>> 16, y = n[(ct | 4) >> 2] | 0, a = y & 65535, 
                y = y >>> 16, p = n[(ct | 8) >> 2] | 0, h = p & 65535, p = p >>> 16, v = n[(ct | 12) >> 2] | 0, 
                o = v & 65535, v = v >>> 16;
                H = n[(ct | 16) >> 2] | 0, S = H & 65535, H = H >>> 16, I = n[(ct | 20) >> 2] | 0, 
                A = I & 65535, I = I >>> 16, q = n[(ct | 24) >> 2] | 0, E = q & 65535, q = q >>> 16, 
                M = n[(ct | 28) >> 2] | 0, _ = M & 65535, M = M >>> 16;
                $ = s(i, S) | 0;
                tt = s(w, S) | 0;
                et = ((s(i, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(w, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                T = et << 16 | $ & 65535;
                $ = (s(i, A) | 0) + (nt & 65535) | 0;
                tt = (s(w, A) | 0) + (nt >>> 16) | 0;
                et = ((s(i, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(w, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                R = et << 16 | $ & 65535;
                $ = (s(i, E) | 0) + (nt & 65535) | 0;
                tt = (s(w, E) | 0) + (nt >>> 16) | 0;
                et = ((s(i, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(w, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                B = et << 16 | $ & 65535;
                $ = (s(i, _) | 0) + (nt & 65535) | 0;
                tt = (s(w, _) | 0) + (nt >>> 16) | 0;
                et = ((s(i, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(w, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                P = et << 16 | $ & 65535;
                K = nt;
                $ = (s(a, S) | 0) + (R & 65535) | 0;
                tt = (s(y, S) | 0) + (R >>> 16) | 0;
                et = ((s(a, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(y, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                R = et << 16 | $ & 65535;
                $ = ((s(a, A) | 0) + (B & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(y, A) | 0) + (B >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(a, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(y, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                B = et << 16 | $ & 65535;
                $ = ((s(a, E) | 0) + (P & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(y, E) | 0) + (P >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(a, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(y, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                P = et << 16 | $ & 65535;
                $ = ((s(a, _) | 0) + (K & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(y, _) | 0) + (K >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(a, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(y, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                K = et << 16 | $ & 65535;
                N = nt;
                $ = (s(h, S) | 0) + (B & 65535) | 0;
                tt = (s(p, S) | 0) + (B >>> 16) | 0;
                et = ((s(h, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(p, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                B = et << 16 | $ & 65535;
                $ = ((s(h, A) | 0) + (P & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(p, A) | 0) + (P >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(h, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(p, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                P = et << 16 | $ & 65535;
                $ = ((s(h, E) | 0) + (K & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(p, E) | 0) + (K >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(h, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(p, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                K = et << 16 | $ & 65535;
                $ = ((s(h, _) | 0) + (N & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(p, _) | 0) + (N >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(h, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(p, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                N = et << 16 | $ & 65535;
                j = nt;
                $ = (s(o, S) | 0) + (P & 65535) | 0;
                tt = (s(v, S) | 0) + (P >>> 16) | 0;
                et = ((s(o, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(v, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                P = et << 16 | $ & 65535;
                $ = ((s(o, A) | 0) + (K & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(v, A) | 0) + (K >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(o, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(v, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                K = et << 16 | $ & 65535;
                $ = ((s(o, E) | 0) + (N & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(v, E) | 0) + (N >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(o, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(v, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                N = et << 16 | $ & 65535;
                $ = ((s(o, _) | 0) + (j & 65535) | 0) + (nt & 65535) | 0;
                tt = ((s(v, _) | 0) + (j >>> 16) | 0) + (nt >>> 16) | 0;
                et = ((s(o, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                nt = ((s(v, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                j = et << 16 | $ & 65535;
                D = nt;
                st = n[(yt | 16) >> 2] | 0;
                $ = (st & 65535) + ((T & 65535) << 1) | 0;
                et = ((st >>> 16) + (T >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 16) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[(yt | 20) >> 2] | 0;
                $ = ((st & 65535) + ((R & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (R >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 20) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[(yt | 24) >> 2] | 0;
                $ = ((st & 65535) + ((B & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (B >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 24) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[(yt | 28) >> 2] | 0;
                $ = ((st & 65535) + ((P & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (P >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[(yt | 28) >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[yt + 32 >> 2] | 0;
                $ = ((st & 65535) + ((K & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (K >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[yt + 32 >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[yt + 36 >> 2] | 0;
                $ = ((st & 65535) + ((N & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (N >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[yt + 36 >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[yt + 40 >> 2] | 0;
                $ = ((st & 65535) + ((j & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (j >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[yt + 40 >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                st = n[yt + 44 >> 2] | 0;
                $ = ((st & 65535) + ((D & 65535) << 1) | 0) + rt | 0;
                et = ((st >>> 16) + (D >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                n[yt + 44 >> 2] = et << 16 | $ & 65535;
                rt = et >>> 16;
                for (lt = 48; !!rt & (lt | 0) < 64; lt = lt + 4 | 0) {
                    st = n[yt + lt >> 2] | 0;
                    $ = (st & 65535) + rt | 0;
                    et = (st >>> 16) + ($ >>> 16) | 0;
                    n[yt + lt >> 2] = et << 16 | $ & 65535;
                    rt = et >>> 16;
                }
            }
            for (at = 32; (at | 0) < (e | 0); at = at << 1) {
                ht = at << 1;
                for (ot = 0; (ot | 0) < (e | 0); ot = ot + ht | 0) {
                    yt = r + (ot << 1) | 0;
                    it = 0;
                    for (ft = 0; (ft | 0) < (at | 0); ft = ft + 32 | 0) {
                        ct = (t + ot | 0) + ft | 0;
                        w = n[ct >> 2] | 0, i = w & 65535, w = w >>> 16, y = n[(ct | 4) >> 2] | 0, a = y & 65535, 
                        y = y >>> 16, p = n[(ct | 8) >> 2] | 0, h = p & 65535, p = p >>> 16, v = n[(ct | 12) >> 2] | 0, 
                        o = v & 65535, v = v >>> 16, d = n[(ct | 16) >> 2] | 0, f = d & 65535, d = d >>> 16, 
                        g = n[(ct | 20) >> 2] | 0, u = g & 65535, g = g >>> 16, b = n[(ct | 24) >> 2] | 0, 
                        l = b & 65535, b = b >>> 16, m = n[(ct | 28) >> 2] | 0, c = m & 65535, m = m >>> 16;
                        G = F = W = V = J = Q = X = Y = rt = 0;
                        for (ut = 0; (ut | 0) < (at | 0); ut = ut + 32 | 0) {
                            wt = ((t + ot | 0) + at | 0) + ut | 0;
                            H = n[wt >> 2] | 0, S = H & 65535, H = H >>> 16, I = n[(wt | 4) >> 2] | 0, A = I & 65535, 
                            I = I >>> 16, q = n[(wt | 8) >> 2] | 0, E = q & 65535, q = q >>> 16, M = n[(wt | 12) >> 2] | 0, 
                            _ = M & 65535, M = M >>> 16, C = n[(wt | 16) >> 2] | 0, k = C & 65535, C = C >>> 16, 
                            Z = n[(wt | 20) >> 2] | 0, U = Z & 65535, Z = Z >>> 16, z = n[(wt | 24) >> 2] | 0, 
                            x = z & 65535, z = z >>> 16, O = n[(wt | 28) >> 2] | 0, L = O & 65535, O = O >>> 16;
                            T = R = B = P = K = N = j = D = 0;
                            $ = ((s(i, S) | 0) + (T & 65535) | 0) + (G & 65535) | 0;
                            tt = ((s(w, S) | 0) + (T >>> 16) | 0) + (G >>> 16) | 0;
                            et = ((s(i, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(w, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            T = et << 16 | $ & 65535;
                            $ = ((s(i, A) | 0) + (R & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(w, A) | 0) + (R >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(i, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(w, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            R = et << 16 | $ & 65535;
                            $ = ((s(i, E) | 0) + (B & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(w, E) | 0) + (B >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(i, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(w, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            B = et << 16 | $ & 65535;
                            $ = ((s(i, _) | 0) + (P & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(w, _) | 0) + (P >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(i, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(w, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            P = et << 16 | $ & 65535;
                            $ = ((s(i, k) | 0) + (K & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(w, k) | 0) + (K >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(i, C) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(w, C) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            K = et << 16 | $ & 65535;
                            $ = ((s(i, U) | 0) + (N & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(w, U) | 0) + (N >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(i, Z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(w, Z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            N = et << 16 | $ & 65535;
                            $ = ((s(i, x) | 0) + (j & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(w, x) | 0) + (j >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(i, z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(w, z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            j = et << 16 | $ & 65535;
                            $ = ((s(i, L) | 0) + (D & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(w, L) | 0) + (D >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(i, O) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(w, O) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            D = et << 16 | $ & 65535;
                            G = nt;
                            $ = ((s(a, S) | 0) + (R & 65535) | 0) + (F & 65535) | 0;
                            tt = ((s(y, S) | 0) + (R >>> 16) | 0) + (F >>> 16) | 0;
                            et = ((s(a, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(y, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            R = et << 16 | $ & 65535;
                            $ = ((s(a, A) | 0) + (B & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(y, A) | 0) + (B >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(a, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(y, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            B = et << 16 | $ & 65535;
                            $ = ((s(a, E) | 0) + (P & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(y, E) | 0) + (P >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(a, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(y, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            P = et << 16 | $ & 65535;
                            $ = ((s(a, _) | 0) + (K & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(y, _) | 0) + (K >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(a, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(y, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            K = et << 16 | $ & 65535;
                            $ = ((s(a, k) | 0) + (N & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(y, k) | 0) + (N >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(a, C) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(y, C) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            N = et << 16 | $ & 65535;
                            $ = ((s(a, U) | 0) + (j & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(y, U) | 0) + (j >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(a, Z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(y, Z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            j = et << 16 | $ & 65535;
                            $ = ((s(a, x) | 0) + (D & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(y, x) | 0) + (D >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(a, z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(y, z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            D = et << 16 | $ & 65535;
                            $ = ((s(a, L) | 0) + (G & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(y, L) | 0) + (G >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(a, O) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(y, O) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            G = et << 16 | $ & 65535;
                            F = nt;
                            $ = ((s(h, S) | 0) + (B & 65535) | 0) + (W & 65535) | 0;
                            tt = ((s(p, S) | 0) + (B >>> 16) | 0) + (W >>> 16) | 0;
                            et = ((s(h, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(p, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            B = et << 16 | $ & 65535;
                            $ = ((s(h, A) | 0) + (P & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(p, A) | 0) + (P >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(h, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(p, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            P = et << 16 | $ & 65535;
                            $ = ((s(h, E) | 0) + (K & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(p, E) | 0) + (K >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(h, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(p, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            K = et << 16 | $ & 65535;
                            $ = ((s(h, _) | 0) + (N & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(p, _) | 0) + (N >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(h, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(p, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            N = et << 16 | $ & 65535;
                            $ = ((s(h, k) | 0) + (j & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(p, k) | 0) + (j >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(h, C) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(p, C) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            j = et << 16 | $ & 65535;
                            $ = ((s(h, U) | 0) + (D & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(p, U) | 0) + (D >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(h, Z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(p, Z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            D = et << 16 | $ & 65535;
                            $ = ((s(h, x) | 0) + (G & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(p, x) | 0) + (G >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(h, z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(p, z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            G = et << 16 | $ & 65535;
                            $ = ((s(h, L) | 0) + (F & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(p, L) | 0) + (F >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(h, O) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(p, O) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            F = et << 16 | $ & 65535;
                            W = nt;
                            $ = ((s(o, S) | 0) + (P & 65535) | 0) + (V & 65535) | 0;
                            tt = ((s(v, S) | 0) + (P >>> 16) | 0) + (V >>> 16) | 0;
                            et = ((s(o, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(v, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            P = et << 16 | $ & 65535;
                            $ = ((s(o, A) | 0) + (K & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(v, A) | 0) + (K >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(o, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(v, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            K = et << 16 | $ & 65535;
                            $ = ((s(o, E) | 0) + (N & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(v, E) | 0) + (N >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(o, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(v, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            N = et << 16 | $ & 65535;
                            $ = ((s(o, _) | 0) + (j & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(v, _) | 0) + (j >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(o, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(v, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            j = et << 16 | $ & 65535;
                            $ = ((s(o, k) | 0) + (D & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(v, k) | 0) + (D >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(o, C) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(v, C) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            D = et << 16 | $ & 65535;
                            $ = ((s(o, U) | 0) + (G & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(v, U) | 0) + (G >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(o, Z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(v, Z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            G = et << 16 | $ & 65535;
                            $ = ((s(o, x) | 0) + (F & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(v, x) | 0) + (F >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(o, z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(v, z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            F = et << 16 | $ & 65535;
                            $ = ((s(o, L) | 0) + (W & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(v, L) | 0) + (W >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(o, O) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(v, O) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            W = et << 16 | $ & 65535;
                            V = nt;
                            $ = ((s(f, S) | 0) + (K & 65535) | 0) + (J & 65535) | 0;
                            tt = ((s(d, S) | 0) + (K >>> 16) | 0) + (J >>> 16) | 0;
                            et = ((s(f, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(d, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            K = et << 16 | $ & 65535;
                            $ = ((s(f, A) | 0) + (N & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(d, A) | 0) + (N >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(f, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(d, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            N = et << 16 | $ & 65535;
                            $ = ((s(f, E) | 0) + (j & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(d, E) | 0) + (j >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(f, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(d, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            j = et << 16 | $ & 65535;
                            $ = ((s(f, _) | 0) + (D & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(d, _) | 0) + (D >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(f, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(d, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            D = et << 16 | $ & 65535;
                            $ = ((s(f, k) | 0) + (G & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(d, k) | 0) + (G >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(f, C) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(d, C) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            G = et << 16 | $ & 65535;
                            $ = ((s(f, U) | 0) + (F & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(d, U) | 0) + (F >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(f, Z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(d, Z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            F = et << 16 | $ & 65535;
                            $ = ((s(f, x) | 0) + (W & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(d, x) | 0) + (W >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(f, z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(d, z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            W = et << 16 | $ & 65535;
                            $ = ((s(f, L) | 0) + (V & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(d, L) | 0) + (V >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(f, O) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(d, O) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            V = et << 16 | $ & 65535;
                            J = nt;
                            $ = ((s(u, S) | 0) + (N & 65535) | 0) + (Q & 65535) | 0;
                            tt = ((s(g, S) | 0) + (N >>> 16) | 0) + (Q >>> 16) | 0;
                            et = ((s(u, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(g, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            N = et << 16 | $ & 65535;
                            $ = ((s(u, A) | 0) + (j & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(g, A) | 0) + (j >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(u, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(g, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            j = et << 16 | $ & 65535;
                            $ = ((s(u, E) | 0) + (D & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(g, E) | 0) + (D >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(u, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(g, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            D = et << 16 | $ & 65535;
                            $ = ((s(u, _) | 0) + (G & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(g, _) | 0) + (G >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(u, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(g, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            G = et << 16 | $ & 65535;
                            $ = ((s(u, k) | 0) + (F & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(g, k) | 0) + (F >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(u, C) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(g, C) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            F = et << 16 | $ & 65535;
                            $ = ((s(u, U) | 0) + (W & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(g, U) | 0) + (W >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(u, Z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(g, Z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            W = et << 16 | $ & 65535;
                            $ = ((s(u, x) | 0) + (V & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(g, x) | 0) + (V >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(u, z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(g, z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            V = et << 16 | $ & 65535;
                            $ = ((s(u, L) | 0) + (J & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(g, L) | 0) + (J >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(u, O) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(g, O) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            J = et << 16 | $ & 65535;
                            Q = nt;
                            $ = ((s(l, S) | 0) + (j & 65535) | 0) + (X & 65535) | 0;
                            tt = ((s(b, S) | 0) + (j >>> 16) | 0) + (X >>> 16) | 0;
                            et = ((s(l, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(b, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            j = et << 16 | $ & 65535;
                            $ = ((s(l, A) | 0) + (D & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(b, A) | 0) + (D >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(l, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(b, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            D = et << 16 | $ & 65535;
                            $ = ((s(l, E) | 0) + (G & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(b, E) | 0) + (G >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(l, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(b, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            G = et << 16 | $ & 65535;
                            $ = ((s(l, _) | 0) + (F & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(b, _) | 0) + (F >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(l, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(b, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            F = et << 16 | $ & 65535;
                            $ = ((s(l, k) | 0) + (W & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(b, k) | 0) + (W >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(l, C) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(b, C) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            W = et << 16 | $ & 65535;
                            $ = ((s(l, U) | 0) + (V & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(b, U) | 0) + (V >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(l, Z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(b, Z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            V = et << 16 | $ & 65535;
                            $ = ((s(l, x) | 0) + (J & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(b, x) | 0) + (J >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(l, z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(b, z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            J = et << 16 | $ & 65535;
                            $ = ((s(l, L) | 0) + (Q & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(b, L) | 0) + (Q >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(l, O) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(b, O) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            Q = et << 16 | $ & 65535;
                            X = nt;
                            $ = ((s(c, S) | 0) + (D & 65535) | 0) + (Y & 65535) | 0;
                            tt = ((s(m, S) | 0) + (D >>> 16) | 0) + (Y >>> 16) | 0;
                            et = ((s(c, H) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(m, H) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            D = et << 16 | $ & 65535;
                            $ = ((s(c, A) | 0) + (G & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(m, A) | 0) + (G >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(c, I) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(m, I) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            G = et << 16 | $ & 65535;
                            $ = ((s(c, E) | 0) + (F & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(m, E) | 0) + (F >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(c, q) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(m, q) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            F = et << 16 | $ & 65535;
                            $ = ((s(c, _) | 0) + (W & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(m, _) | 0) + (W >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(c, M) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(m, M) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            W = et << 16 | $ & 65535;
                            $ = ((s(c, k) | 0) + (V & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(m, k) | 0) + (V >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(c, C) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(m, C) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            V = et << 16 | $ & 65535;
                            $ = ((s(c, U) | 0) + (J & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(m, U) | 0) + (J >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(c, Z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(m, Z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            J = et << 16 | $ & 65535;
                            $ = ((s(c, x) | 0) + (Q & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(m, x) | 0) + (Q >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(c, z) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(m, z) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            Q = et << 16 | $ & 65535;
                            $ = ((s(c, L) | 0) + (X & 65535) | 0) + (nt & 65535) | 0;
                            tt = ((s(m, L) | 0) + (X >>> 16) | 0) + (nt >>> 16) | 0;
                            et = ((s(c, O) | 0) + (tt & 65535) | 0) + ($ >>> 16) | 0;
                            nt = ((s(m, O) | 0) + (tt >>> 16) | 0) + (et >>> 16) | 0;
                            X = et << 16 | $ & 65535;
                            Y = nt;
                            lt = at + (ft + ut | 0) | 0;
                            st = n[yt + lt >> 2] | 0;
                            $ = ((st & 65535) + ((T & 65535) << 1) | 0) + rt | 0;
                            et = ((st >>> 16) + (T >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                            n[yt + lt >> 2] = et << 16 | $ & 65535;
                            rt = et >>> 16;
                            lt = lt + 4 | 0;
                            st = n[yt + lt >> 2] | 0;
                            $ = ((st & 65535) + ((R & 65535) << 1) | 0) + rt | 0;
                            et = ((st >>> 16) + (R >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                            n[yt + lt >> 2] = et << 16 | $ & 65535;
                            rt = et >>> 16;
                            lt = lt + 4 | 0;
                            st = n[yt + lt >> 2] | 0;
                            $ = ((st & 65535) + ((B & 65535) << 1) | 0) + rt | 0;
                            et = ((st >>> 16) + (B >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                            n[yt + lt >> 2] = et << 16 | $ & 65535;
                            rt = et >>> 16;
                            lt = lt + 4 | 0;
                            st = n[yt + lt >> 2] | 0;
                            $ = ((st & 65535) + ((P & 65535) << 1) | 0) + rt | 0;
                            et = ((st >>> 16) + (P >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                            n[yt + lt >> 2] = et << 16 | $ & 65535;
                            rt = et >>> 16;
                            lt = lt + 4 | 0;
                            st = n[yt + lt >> 2] | 0;
                            $ = ((st & 65535) + ((K & 65535) << 1) | 0) + rt | 0;
                            et = ((st >>> 16) + (K >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                            n[yt + lt >> 2] = et << 16 | $ & 65535;
                            rt = et >>> 16;
                            lt = lt + 4 | 0;
                            st = n[yt + lt >> 2] | 0;
                            $ = ((st & 65535) + ((N & 65535) << 1) | 0) + rt | 0;
                            et = ((st >>> 16) + (N >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                            n[yt + lt >> 2] = et << 16 | $ & 65535;
                            rt = et >>> 16;
                            lt = lt + 4 | 0;
                            st = n[yt + lt >> 2] | 0;
                            $ = ((st & 65535) + ((j & 65535) << 1) | 0) + rt | 0;
                            et = ((st >>> 16) + (j >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                            n[yt + lt >> 2] = et << 16 | $ & 65535;
                            rt = et >>> 16;
                            lt = lt + 4 | 0;
                            st = n[yt + lt >> 2] | 0;
                            $ = ((st & 65535) + ((D & 65535) << 1) | 0) + rt | 0;
                            et = ((st >>> 16) + (D >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                            n[yt + lt >> 2] = et << 16 | $ & 65535;
                            rt = et >>> 16;
                        }
                        lt = at + (ft + ut | 0) | 0;
                        st = n[yt + lt >> 2] | 0;
                        $ = (((st & 65535) + ((G & 65535) << 1) | 0) + rt | 0) + it | 0;
                        et = ((st >>> 16) + (G >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                        n[yt + lt >> 2] = et << 16 | $ & 65535;
                        rt = et >>> 16;
                        lt = lt + 4 | 0;
                        st = n[yt + lt >> 2] | 0;
                        $ = ((st & 65535) + ((F & 65535) << 1) | 0) + rt | 0;
                        et = ((st >>> 16) + (F >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                        n[yt + lt >> 2] = et << 16 | $ & 65535;
                        rt = et >>> 16;
                        lt = lt + 4 | 0;
                        st = n[yt + lt >> 2] | 0;
                        $ = ((st & 65535) + ((W & 65535) << 1) | 0) + rt | 0;
                        et = ((st >>> 16) + (W >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                        n[yt + lt >> 2] = et << 16 | $ & 65535;
                        rt = et >>> 16;
                        lt = lt + 4 | 0;
                        st = n[yt + lt >> 2] | 0;
                        $ = ((st & 65535) + ((V & 65535) << 1) | 0) + rt | 0;
                        et = ((st >>> 16) + (V >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                        n[yt + lt >> 2] = et << 16 | $ & 65535;
                        rt = et >>> 16;
                        lt = lt + 4 | 0;
                        st = n[yt + lt >> 2] | 0;
                        $ = ((st & 65535) + ((J & 65535) << 1) | 0) + rt | 0;
                        et = ((st >>> 16) + (J >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                        n[yt + lt >> 2] = et << 16 | $ & 65535;
                        rt = et >>> 16;
                        lt = lt + 4 | 0;
                        st = n[yt + lt >> 2] | 0;
                        $ = ((st & 65535) + ((Q & 65535) << 1) | 0) + rt | 0;
                        et = ((st >>> 16) + (Q >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                        n[yt + lt >> 2] = et << 16 | $ & 65535;
                        rt = et >>> 16;
                        lt = lt + 4 | 0;
                        st = n[yt + lt >> 2] | 0;
                        $ = ((st & 65535) + ((X & 65535) << 1) | 0) + rt | 0;
                        et = ((st >>> 16) + (X >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                        n[yt + lt >> 2] = et << 16 | $ & 65535;
                        rt = et >>> 16;
                        lt = lt + 4 | 0;
                        st = n[yt + lt >> 2] | 0;
                        $ = ((st & 65535) + ((Y & 65535) << 1) | 0) + rt | 0;
                        et = ((st >>> 16) + (Y >>> 16 << 1) | 0) + ($ >>> 16) | 0;
                        n[yt + lt >> 2] = et << 16 | $ & 65535;
                        it = et >>> 16;
                    }
                    for (lt = lt + 4 | 0; !!it & (lt | 0) < ht << 1; lt = lt + 4 | 0) {
                        st = n[yt + lt >> 2] | 0;
                        $ = (st & 65535) + it | 0;
                        et = (st >>> 16) + ($ >>> 16) | 0;
                        n[yt + lt >> 2] = et << 16 | $ & 65535;
                        it = et >>> 16;
                    }
                }
            }
        }
        function g(t, e, r, i, a, h) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            a = a | 0;
            h = h | 0;
            var o = 0, u = 0, l = 0, c = 0, w = 0, y = 0, p = 0, v = 0, d = 0, g = 0, b = 0, m = 0, S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0;
            f(e, t, a);
            for (k = e - 1 & -4; (k | 0) >= 0; k = k - 4 | 0) {
                o = n[t + k >> 2] | 0;
                if (o) {
                    e = k;
                    break;
                }
            }
            for (k = i - 1 & -4; (k | 0) >= 0; k = k - 4 | 0) {
                u = n[r + k >> 2] | 0;
                if (u) {
                    i = k;
                    break;
                }
            }
            while ((u & 2147483648) == 0) {
                u = u << 1;
                l = l + 1 | 0;
            }
            w = n[t + e >> 2] | 0;
            if (l) c = w >>> (32 - l | 0);
            for (k = e - 4 | 0; (k | 0) >= 0; k = k - 4 | 0) {
                o = n[t + k >> 2] | 0;
                n[a + k + 4 >> 2] = w << l | (l ? o >>> (32 - l | 0) : 0);
                w = o;
            }
            n[a >> 2] = w << l;
            if (l) {
                y = n[r + i >> 2] | 0;
                for (k = i - 4 | 0; (k | 0) >= 0; k = k - 4 | 0) {
                    u = n[r + k >> 2] | 0;
                    n[r + k + 4 >> 2] = y << l | u >>> (32 - l | 0);
                    y = u;
                }
                n[r >> 2] = y << l;
            }
            y = n[r + i >> 2] | 0;
            p = y >>> 16, v = y & 65535;
            for (k = e; (k | 0) >= (i | 0); k = k - 4 | 0) {
                U = k - i | 0;
                w = n[a + k >> 2] | 0;
                d = (c >>> 0) / (p >>> 0) | 0, b = (c >>> 0) % (p >>> 0) | 0, S = s(d, v) | 0;
                while ((d | 0) == 65536 | S >>> 0 > (b << 16 | w >>> 16) >>> 0) {
                    d = d - 1 | 0, b = b + p | 0, S = S - v | 0;
                    if ((b | 0) >= 65536) break;
                }
                E = 0, _ = 0;
                for (x = 0; (x | 0) <= (i | 0); x = x + 4 | 0) {
                    u = n[r + x >> 2] | 0;
                    S = (s(d, u & 65535) | 0) + (E >>> 16) | 0;
                    A = (s(d, u >>> 16) | 0) + (S >>> 16) | 0;
                    u = E & 65535 | S << 16;
                    E = A;
                    o = n[a + U + x >> 2] | 0;
                    S = ((o & 65535) - (u & 65535) | 0) + _ | 0;
                    A = ((o >>> 16) - (u >>> 16) | 0) + (S >> 16) | 0;
                    n[a + U + x >> 2] = A << 16 | S & 65535;
                    _ = A >> 16;
                }
                S = ((c & 65535) - (E & 65535) | 0) + _ | 0;
                A = ((c >>> 16) - (E >>> 16) | 0) + (S >> 16) | 0;
                n[a + U + x >> 2] = c = A << 16 | S & 65535;
                _ = A >> 16;
                if (_) {
                    d = d - 1 | 0, b = b - p | 0;
                    _ = 0;
                    for (x = 0; (x | 0) <= (i | 0); x = x + 4 | 0) {
                        u = n[r + x >> 2] | 0;
                        o = n[a + U + x >> 2] | 0;
                        S = ((o & 65535) + (u & 65535) | 0) + _ | 0;
                        A = ((o >>> 16) + (u >>> 16) | 0) + (S >>> 16) | 0;
                        n[a + U + x >> 2] = A << 16 | S & 65535;
                        _ = A >>> 16;
                    }
                    n[a + U + x >> 2] = c = c + _ | 0;
                }
                w = n[a + k >> 2] | 0;
                o = c << 16 | w >>> 16;
                g = (o >>> 0) / (p >>> 0) | 0, m = (o >>> 0) % (p >>> 0) | 0, S = s(g, v) | 0;
                while ((g | 0) == 65536 | S >>> 0 > (m << 16 | w & 65535) >>> 0) {
                    g = g - 1 | 0, m = m + p | 0, S = S - v | 0;
                    if ((m | 0) >= 65536) break;
                }
                E = 0, _ = 0;
                for (x = 0; (x | 0) <= (i | 0); x = x + 4 | 0) {
                    u = n[r + x >> 2] | 0;
                    S = (s(g, u & 65535) | 0) + (E & 65535) | 0;
                    A = ((s(g, u >>> 16) | 0) + (S >>> 16) | 0) + (E >>> 16) | 0;
                    u = S & 65535 | A << 16;
                    E = A >>> 16;
                    o = n[a + U + x >> 2] | 0;
                    S = ((o & 65535) - (u & 65535) | 0) + _ | 0;
                    A = ((o >>> 16) - (u >>> 16) | 0) + (S >> 16) | 0;
                    _ = A >> 16;
                    n[a + U + x >> 2] = A << 16 | S & 65535;
                }
                S = ((c & 65535) - (E & 65535) | 0) + _ | 0;
                A = ((c >>> 16) - (E >>> 16) | 0) + (S >> 16) | 0;
                n[a + U + x >> 2] = c = A << 16 | S & 65535;
                _ = A >> 16;
                if (_) {
                    g = g - 1 | 0, m = m + p | 0;
                    _ = 0;
                    for (x = 0; (x | 0) <= (i | 0); x = x + 4 | 0) {
                        u = n[r + x >> 2] | 0;
                        o = n[a + U + x >> 2] | 0;
                        S = ((o & 65535) + (u & 65535) | 0) + _ | 0;
                        A = ((o >>> 16) + (u >>> 16) | 0) + (S >>> 16) | 0;
                        _ = A >>> 16;
                        n[a + U + x >> 2] = S & 65535 | A << 16;
                    }
                    n[a + U + x >> 2] = c + _ | 0;
                }
                n[h + U >> 2] = d << 16 | g;
                c = n[a + k >> 2] | 0;
            }
            if (l) {
                w = n[a >> 2] | 0;
                for (k = 4; (k | 0) <= (i | 0); k = k + 4 | 0) {
                    o = n[a + k >> 2] | 0;
                    n[a + k - 4 >> 2] = o << (32 - l | 0) | w >>> l;
                    w = o;
                }
                n[a + i >> 2] = w >>> l;
            }
        }
        function b(t, e, r, i, a, l) {
            t = t | 0;
            e = e | 0;
            r = r | 0;
            i = i | 0;
            a = a | 0;
            l = l | 0;
            var w = 0, y = 0, v = 0, d = 0, g = 0, b = 0, m = 0, S = 0, A = 0, E = 0, _ = 0, k = 0, U = 0, x = 0;
            w = h(i << 1) | 0;
            u(i << 1, 0, w);
            f(e, t, w);
            for (k = 0; (k | 0) < (i | 0); k = k + 4 | 0) {
                v = n[w + k >> 2] | 0, d = v & 65535, v = v >>> 16;
                b = a >>> 16, g = a & 65535;
                m = s(d, g) | 0, S = ((s(d, b) | 0) + (s(v, g) | 0) | 0) + (m >>> 16) | 0;
                d = m & 65535, v = S & 65535;
                _ = 0;
                for (U = 0; (U | 0) < (i | 0); U = U + 4 | 0) {
                    x = k + U | 0;
                    b = n[r + U >> 2] | 0, g = b & 65535, b = b >>> 16;
                    E = n[w + x >> 2] | 0;
                    m = ((s(d, g) | 0) + (_ & 65535) | 0) + (E & 65535) | 0;
                    S = ((s(d, b) | 0) + (_ >>> 16) | 0) + (E >>> 16) | 0;
                    A = ((s(v, g) | 0) + (S & 65535) | 0) + (m >>> 16) | 0;
                    _ = ((s(v, b) | 0) + (A >>> 16) | 0) + (S >>> 16) | 0;
                    E = A << 16 | m & 65535;
                    n[w + x >> 2] = E;
                }
                x = k + U | 0;
                E = n[w + x >> 2] | 0;
                m = ((E & 65535) + (_ & 65535) | 0) + y | 0;
                S = ((E >>> 16) + (_ >>> 16) | 0) + (m >>> 16) | 0;
                n[w + x >> 2] = S << 16 | m & 65535;
                y = S >>> 16;
            }
            f(i, w + i | 0, l);
            o(i << 1);
            if (y | (c(r, i, l, i) | 0) <= 0) {
                p(l, i, r, i, l, i) | 0;
            }
        }
        return {
            sreset: a,
            salloc: h,
            sfree: o,
            z: u,
            tst: w,
            neg: l,
            cmp: c,
            add: y,
            sub: p,
            mul: v,
            sqr: d,
            div: g,
            mredc: b
        };
    }
    function ot(t) {
        return t instanceof ft;
    }
    function ft(t) {
        var e = Ae, r = 0, i = 0;
        if (l(t) && (t = s(t)), c(t) && (t = new Uint8Array(t)), void 0 === t) ; else if (u(t)) {
            var n = Math.abs(t);
            n > 4294967295 ? ((e = new Uint32Array(2))[0] = 0 | n, e[1] = n / 4294967296 | 0, 
            r = 52) : n > 0 ? ((e = new Uint32Array(1))[0] = n, r = 32) : (e = Ae, r = 0), i = 0 > t ? -1 : 1;
        } else if (w(t)) {
            if (!(r = 8 * t.length)) return _e;
            e = new Uint32Array(r + 31 >> 5);
            for (var a = t.length - 4; a >= 0; a -= 4) e[t.length - 4 - a >> 2] = t[a] << 24 | t[a + 1] << 16 | t[a + 2] << 8 | t[a + 3];
            -3 === a ? e[e.length - 1] = t[0] : -2 === a ? e[e.length - 1] = t[0] << 8 | t[1] : -1 === a && (e[e.length - 1] = t[0] << 16 | t[1] << 8 | t[2]), 
            i = 1;
        } else {
            if ("object" != typeof t || null === t) throw new TypeError("number is of unexpected type");
            e = new Uint32Array(t.limbs), r = t.bitLength, i = t.sign;
        }
        this.limbs = e, this.bitLength = r, this.sign = i;
    }
    function ut(t, e) {
        ot(t) || (t = new ft(t)), ot(e) || (e = new ft(e));
        var r = t.sign, i = e.sign;
        0 > r && (t = t.negate()), 0 > i && (e = e.negate());
        var n = t.compare(e);
        if (0 > n) {
            var s = t;
            t = e, e = s, s = r, r = i, i = s;
        }
        var a, h, o, f = ke, u = _e, l = e.bitLength, c = _e, w = ke, y = t.bitLength;
        for (a = t.divide(e); (h = a.remainder) !== _e; ) o = a.quotient, a = f.subtract(o.multiply(u).clamp(l)).clamp(l), 
        f = u, u = a, a = c.subtract(o.multiply(w).clamp(y)).clamp(y), c = w, w = a, t = e, 
        e = h, a = t.divide(e);
        if (0 > r && (u = u.negate()), 0 > i && (w = w.negate()), 0 > n) {
            s = u;
            u = w, w = s;
        }
        return {
            gcd: e,
            x: u,
            y: w
        };
    }
    function lt() {
        if (ft.apply(this, arguments), this.valueOf() < 1) throw new RangeError();
        var t;
        if (!(this.bitLength <= 32) && 1 & this.limbs[0]) {
            var e = 1 + (this.bitLength + 31 & -32), r = new Uint32Array(e + 31 >> 5);
            r[r.length - 1] = 1, (t = new ft()).sign = 1, t.bitLength = e, t.limbs = r;
            var i = function(t, e) {
                var r, i, n, s, a = 0 > t ? -1 : 1, h = 0 > e ? -1 : 1, o = 1, f = 0, u = 0, l = 1;
                for ((s = (e *= h) > (t *= a)) && (n = t, t = e, e = n, n = a, a = h, h = n), r = t - (i = Math.floor(t / e)) * e; r; ) n = o - i * f, 
                o = f, f = n, n = u - i * l, u = l, l = n, t = e, e = r, r = t - (i = Math.floor(t / e)) * e;
                return f *= a, l *= h, s && (n = f, f = l, l = n), {
                    gcd: e,
                    x: f,
                    y: l
                };
            }(4294967296, this.limbs[0]).y;
            this.coefficient = 0 > i ? -i : 4294967296 - i, this.comodulus = t, this.comodulusRemainder = t.divide(this).remainder, 
            this.comodulusRemainderSquare = t.square().divide(this).remainder;
        }
    }
    function ct(t, e) {
        var r = t.limbs, i = r.length, n = e.limbs, s = n.length, a = e.coefficient;
        ht.sreset();
        var h = ht.salloc(i << 2), o = ht.salloc(s << 2), f = ht.salloc(s << 2);
        ht.z(f - h + (s << 2), 0, h), Se.set(r, h >> 2), Se.set(n, o >> 2), ht.mredc(h, i << 2, o, s << 2, a, f);
        var u = new ft();
        return u.limbs = new Uint32Array(Se.subarray(f >> 2, (f >> 2) + s)), u.bitLength = e.bitLength, 
        u.sign = 1, u;
    }
    function wt(t) {
        var e = new ft(this), r = 0;
        for (e.limbs[0] -= 1; 0 === e.limbs[r >> 5]; ) r += 32;
        for (;0 == (e.limbs[r >> 5] >> (31 & r) & 1); ) r++;
        e = e.slice(r);
        for (var i = new lt(this), n = this.subtract(ke), s = new ft(this), a = this.limbs.length - 1; 0 === s.limbs[a]; ) a--;
        for (;--t >= 0; ) {
            for (nt(s.limbs), s.limbs[0] < 2 && (s.limbs[0] += 2); s.compare(n) >= 0; ) s.limbs[a] >>>= 1;
            var h = i.power(s, e);
            if (0 !== h.compare(ke) && 0 !== h.compare(n)) {
                for (var o = r; --o > 0; ) {
                    if (0 === (h = h.square().divide(i).remainder).compare(ke)) return !1;
                    if (0 === h.compare(n)) break;
                }
                if (0 === o) return !1;
            }
        }
        return !0;
    }
    function yt(t, r) {
        var i = t + 31 >> 5, n = new ft({
            sign: 1,
            bitLength: t,
            limbs: i
        }), s = n.limbs, a = 1e4;
        512 >= t && (a = 2200), 256 >= t && (a = 600);
        var h = function(t) {
            if (xe.length >= t) return xe.slice(0, t);
            for (var e = xe[xe.length - 1] + 2; xe.length < t; e += 2) {
                for (var r = 0, i = xe[r]; e >= i * i && e % i != 0; i = xe[++r]) ;
                i * i > e && xe.push(e);
            }
            return xe;
        }(a), o = new Uint32Array(a), u = t * e.Math.LN2 | 0, l = 27;
        for (t >= 250 && (l = 12), t >= 450 && (l = 6), t >= 850 && (l = 3), t >= 1300 && (l = 2); ;) {
            nt(s), s[0] |= 1, s[i - 1] |= 1 << (t - 1 & 31), 31 & t && (s[i - 1] &= f(t + 1 & 31) - 1), 
            o[0] = 1;
            for (var c = 1; a > c; c++) o[c] = n.divide(h[c]).remainder.valueOf();
            t: for (var w = 0; u > w; w += 2, s[0] += 2) {
                for (c = 1; a > c; c++) if ((o[c] + w) % h[c] == 0) continue t;
                if (("function" != typeof r || r(n)) && wt.call(n, l)) return n;
            }
        }
    }
    function pt(t) {
        t = t || {}, this.key = null, this.result = null, this.reset(t);
    }
    function vt(t) {
        t = t || {}, this.result = null;
        var e = t.key;
        if (void 0 !== e) {
            if (!(e instanceof Array)) throw new TypeError("unexpected key type");
            var r = e.length;
            if (2 !== r && 3 !== r && 8 !== r) throw new SyntaxError("unexpected key type");
            var i = [];
            i[0] = new lt(e[0]), i[1] = new ft(e[1]), r > 2 && (i[2] = new ft(e[2])), r > 3 && (i[3] = new lt(e[3]), 
            i[4] = new lt(e[4]), i[5] = new ft(e[5]), i[6] = new ft(e[6]), i[7] = new ft(e[7])), 
            this.key = i;
        }
        return this;
    }
    function dt(t) {
        if (!this.key) throw new r("no key is associated with the instance");
        var e;
        if (l(t) && (t = s(t)), c(t) && (t = new Uint8Array(t)), w(t)) e = new ft(t); else {
            if (!ot(t)) throw new TypeError("unexpected data type");
            e = t;
        }
        if (this.key[0].compare(e) <= 0) throw new RangeError("data too large");
        var i = this.key[0], n = this.key[1], a = i.power(e, n).toBytes(), h = i.bitLength + 7 >> 3;
        if (a.length < h) {
            var o = new Uint8Array(h);
            o.set(a, h - a.length), a = o;
        }
        return this.result = a, this;
    }
    function gt(t) {
        if (!this.key) throw new r("no key is associated with the instance");
        if (this.key.length < 3) throw new r("key isn't suitable for decription");
        var e, i;
        if (l(t) && (t = s(t)), c(t) && (t = new Uint8Array(t)), w(t)) e = new ft(t); else {
            if (!ot(t)) throw new TypeError("unexpected data type");
            e = t;
        }
        if (this.key[0].compare(e) <= 0) throw new RangeError("data too large");
        if (this.key.length > 3) {
            for (var n = this.key[0], a = this.key[3], h = this.key[4], o = this.key[5], f = this.key[6], u = this.key[7], y = a.power(e, o), p = h.power(e, f), v = y.subtract(p); v.sign < 0; ) v = v.add(a);
            i = a.reduce(u.multiply(v)).multiply(h).add(p).clamp(n.bitLength).toBytes();
        } else {
            n = this.key[0];
            var d = this.key[2];
            i = n.power(e, d).toBytes();
        }
        var g = n.bitLength + 7 >> 3;
        if (i.length < g) {
            var b = new Uint8Array(g);
            b.set(i, g - i.length), i = b;
        }
        return this.result = i, this;
    }
    function bt(t, e) {
        if (t = t || 2048, e = e || 65537, 512 > t) throw new i("bit length is too small");
        if (l(e) && (e = s(e)), c(e) && (e = new Uint8Array(e)), !(w(e) || u(e) || ot(e))) throw new TypeError("unexpected exponent type");
        if (0 == (1 & (e = new ft(e)).limbs[0])) throw new i("exponent must be an odd number");
        var r, n, a, h, o, f, y, p;
        a = yt(t >> 1, function(t) {
            return (o = new ft(t)).limbs[0] -= 1, 1 == ut(o, e).gcd.valueOf();
        }), h = yt(t - (t >> 1), function(i) {
            return !!((r = new lt(a.multiply(i))).limbs[(t + 31 >> 5) - 1] >>> (t - 1 & 31)) && ((f = new ft(i)).limbs[0] -= 1, 
            1 == ut(f, e).gcd.valueOf());
        }), y = (n = new lt(o.multiply(f)).inverse(e)).divide(o).remainder, p = n.divide(f).remainder, 
        a = new lt(a), h = new lt(h);
        var v = a.inverse(h);
        return [ r, e, n, a, h, y, p, v ];
    }
    function mt(t) {
        if (!(t = t || {}).hash) throw new SyntaxError("option 'hash' is required");
        if (!t.hash.HASH_SIZE) throw new SyntaxError("option 'hash' supplied doesn't seem to be a valid hash function");
        this.hash = t.hash, this.label = null, this.reset(t);
    }
    function St(t, e) {
        t = t || "", e = e || 0;
        for (var r = this.hash.HASH_SIZE, i = new Uint8Array(e), n = new Uint8Array(4), s = Math.ceil(e / r), a = 0; s > a; a++) {
            n[0] = a >>> 24, n[1] = a >>> 16 & 255, n[2] = a >>> 8 & 255, n[3] = 255 & a;
            var h = i.subarray(a * r), o = this.hash.reset().process(t).process(n).finish().result;
            o.length > h.length && (o = o.subarray(0, h.length)), h.set(o);
        }
        return i;
    }
    function At(t) {
        if (!(t = t || {}).hash) throw new SyntaxError("option 'hash' is required");
        if (!t.hash.HASH_SIZE) throw new SyntaxError("option 'hash' supplied doesn't seem to be a valid hash function");
        this.hash = t.hash, this.saltLength = 4, this.reset(t);
    }
    e.asmCrypto = t, r.prototype = Object.create(Error.prototype, {
        name: {
            value: "IllegalStateError"
        }
    }), i.prototype = Object.create(Error.prototype, {
        name: {
            value: "IllegalArgumentError"
        }
    }), n.prototype = Object.create(Error.prototype, {
        name: {
            value: "SecurityError"
        }
    });
    var Et = e.Float64Array || e.Float32Array;
    t.string_to_bytes = s, t.hex_to_bytes = function(t) {
        var e, r = [], i = t.length;
        for (1 & i && (t = "0" + t, i++), e = 0; i > e; e += 2) r.push(parseInt(t.substr(e, 2), 16));
        return new Uint8Array(r);
    }, t.base64_to_bytes = function(t) {
        return s(atob(t));
    }, t.bytes_to_string = a, t.bytes_to_hex = h, t.bytes_to_base64 = o, e.IllegalStateError = r, 
    e.IllegalArgumentError = i, e.SecurityError = n;
    var _t = [ 99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254, 215, 171, 118, 202, 130, 201, 125, 250, 89, 71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 183, 253, 147, 38, 54, 63, 247, 204, 52, 165, 229, 241, 113, 216, 49, 21, 4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39, 178, 117, 9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41, 227, 47, 132, 83, 209, 0, 237, 32, 252, 177, 91, 106, 203, 190, 57, 74, 76, 88, 207, 208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2, 127, 80, 60, 159, 168, 81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16, 255, 243, 210, 205, 12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61, 100, 93, 25, 115, 96, 129, 79, 220, 34, 42, 144, 136, 70, 238, 184, 20, 222, 94, 11, 219, 224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98, 145, 149, 228, 121, 231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234, 101, 122, 174, 8, 186, 120, 37, 46, 28, 166, 180, 198, 232, 221, 116, 31, 75, 189, 139, 138, 112, 62, 181, 102, 72, 3, 246, 14, 97, 53, 87, 185, 134, 193, 29, 158, 225, 248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233, 206, 85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104, 65, 153, 45, 15, 176, 84, 187, 22, 82, 9, 106, 213, 48, 54, 165, 56, 191, 64, 163, 158, 129, 243, 215, 251, 124, 227, 57, 130, 155, 47, 255, 135, 52, 142, 67, 68, 196, 222, 233, 203, 84, 123, 148, 50, 166, 194, 35, 61, 238, 76, 149, 11, 66, 250, 195, 78, 8, 46, 161, 102, 40, 217, 36, 178, 118, 91, 162, 73, 109, 139, 209, 37, 114, 248, 246, 100, 134, 104, 152, 22, 212, 164, 92, 204, 93, 101, 182, 146, 108, 112, 72, 80, 253, 237, 185, 218, 94, 21, 70, 87, 167, 141, 157, 132, 144, 216, 171, 0, 140, 188, 211, 10, 247, 228, 88, 5, 184, 179, 69, 6, 208, 44, 30, 143, 202, 63, 15, 2, 193, 175, 189, 3, 1, 19, 138, 107, 58, 145, 17, 65, 79, 103, 220, 234, 151, 242, 207, 206, 240, 180, 230, 115, 150, 172, 116, 34, 231, 173, 53, 133, 226, 249, 55, 232, 28, 117, 223, 110, 71, 241, 26, 113, 29, 41, 197, 137, 111, 183, 98, 14, 170, 24, 190, 27, 252, 86, 62, 75, 198, 210, 121, 32, 154, 219, 192, 254, 120, 205, 90, 244, 31, 221, 168, 51, 136, 7, 199, 49, 177, 18, 16, 89, 39, 128, 236, 95, 96, 81, 127, 169, 25, 181, 74, 13, 45, 229, 122, 159, 147, 201, 156, 239, 160, 224, 59, 77, 174, 42, 245, 176, 200, 235, 187, 60, 131, 83, 153, 97, 23, 43, 4, 126, 186, 119, 214, 38, 225, 105, 20, 99, 85, 33, 12, 125, 198, 248, 238, 246, 255, 214, 222, 145, 96, 2, 206, 86, 231, 181, 77, 236, 143, 31, 137, 250, 239, 178, 142, 251, 65, 179, 95, 69, 35, 83, 228, 155, 117, 225, 61, 76, 108, 126, 245, 131, 104, 81, 209, 249, 226, 171, 98, 42, 8, 149, 70, 157, 48, 55, 10, 47, 14, 36, 27, 223, 205, 78, 127, 234, 18, 29, 88, 52, 54, 220, 180, 91, 164, 118, 183, 125, 82, 221, 94, 19, 166, 185, 0, 193, 64, 227, 121, 182, 212, 141, 103, 114, 148, 152, 176, 133, 187, 197, 79, 237, 134, 154, 102, 17, 138, 233, 4, 254, 160, 120, 37, 75, 162, 93, 128, 5, 63, 33, 112, 241, 99, 119, 175, 66, 32, 229, 253, 191, 129, 24, 38, 195, 190, 53, 136, 46, 147, 85, 252, 122, 200, 186, 50, 230, 192, 25, 158, 163, 68, 84, 59, 11, 140, 199, 107, 40, 167, 188, 22, 173, 219, 100, 116, 20, 146, 12, 72, 184, 159, 189, 67, 196, 57, 49, 211, 242, 213, 139, 110, 218, 1, 177, 156, 73, 216, 172, 243, 207, 202, 244, 71, 16, 111, 240, 74, 92, 56, 87, 115, 151, 203, 161, 232, 62, 150, 97, 13, 15, 224, 124, 113, 204, 144, 6, 247, 28, 194, 106, 174, 105, 23, 153, 58, 39, 217, 235, 43, 34, 210, 169, 7, 51, 45, 60, 21, 201, 135, 170, 80, 165, 3, 89, 9, 26, 101, 215, 132, 208, 130, 41, 90, 30, 123, 168, 109, 44, 165, 132, 153, 141, 13, 189, 177, 84, 80, 3, 169, 125, 25, 98, 230, 154, 69, 157, 64, 135, 21, 235, 201, 11, 236, 103, 253, 234, 191, 247, 150, 91, 194, 28, 174, 106, 90, 65, 2, 79, 92, 244, 52, 8, 147, 115, 83, 63, 12, 82, 101, 94, 40, 161, 15, 181, 9, 54, 155, 61, 38, 105, 205, 159, 27, 158, 116, 46, 45, 178, 238, 251, 246, 77, 97, 206, 123, 62, 113, 151, 245, 104, 0, 44, 96, 31, 200, 237, 190, 70, 217, 75, 222, 212, 232, 74, 107, 42, 229, 22, 197, 215, 85, 148, 207, 16, 6, 129, 240, 68, 186, 227, 243, 254, 192, 138, 173, 188, 72, 4, 223, 193, 117, 99, 48, 26, 14, 109, 76, 20, 53, 47, 225, 162, 204, 57, 87, 242, 130, 71, 172, 231, 43, 149, 160, 152, 209, 127, 102, 126, 171, 131, 202, 41, 211, 60, 121, 226, 29, 118, 59, 86, 78, 30, 219, 10, 108, 228, 93, 110, 239, 166, 168, 164, 55, 139, 50, 67, 89, 183, 140, 100, 210, 224, 180, 250, 7, 37, 175, 142, 233, 24, 213, 136, 111, 114, 36, 241, 199, 81, 35, 124, 156, 33, 221, 220, 134, 133, 144, 66, 196, 170, 216, 5, 1, 18, 163, 95, 249, 208, 145, 88, 39, 185, 56, 19, 179, 51, 187, 112, 137, 167, 182, 34, 146, 32, 73, 255, 120, 122, 143, 248, 128, 23, 218, 49, 198, 184, 195, 176, 119, 17, 203, 252, 214, 58, 0, 9, 18, 27, 36, 45, 54, 63, 72, 65, 90, 83, 108, 101, 126, 119, 144, 153, 130, 139, 180, 189, 166, 175, 216, 209, 202, 195, 252, 245, 238, 231, 59, 50, 41, 32, 31, 22, 13, 4, 115, 122, 97, 104, 87, 94, 69, 76, 171, 162, 185, 176, 143, 134, 157, 148, 227, 234, 241, 248, 199, 206, 213, 220, 118, 127, 100, 109, 82, 91, 64, 73, 62, 55, 44, 37, 26, 19, 8, 1, 230, 239, 244, 253, 194, 203, 208, 217, 174, 167, 188, 181, 138, 131, 152, 145, 77, 68, 95, 86, 105, 96, 123, 114, 5, 12, 23, 30, 33, 40, 51, 58, 221, 212, 207, 198, 249, 240, 235, 226, 149, 156, 135, 142, 177, 184, 163, 170, 236, 229, 254, 247, 200, 193, 218, 211, 164, 173, 182, 191, 128, 137, 146, 155, 124, 117, 110, 103, 88, 81, 74, 67, 52, 61, 38, 47, 16, 25, 2, 11, 215, 222, 197, 204, 243, 250, 225, 232, 159, 150, 141, 132, 187, 178, 169, 160, 71, 78, 85, 92, 99, 106, 113, 120, 15, 6, 29, 20, 43, 34, 57, 48, 154, 147, 136, 129, 190, 183, 172, 165, 210, 219, 192, 201, 246, 255, 228, 237, 10, 3, 24, 17, 46, 39, 60, 53, 66, 75, 80, 89, 102, 111, 116, 125, 161, 168, 179, 186, 133, 140, 151, 158, 233, 224, 251, 242, 205, 196, 223, 214, 49, 56, 35, 42, 21, 28, 7, 14, 121, 112, 107, 98, 93, 84, 79, 70, 0, 11, 22, 29, 44, 39, 58, 49, 88, 83, 78, 69, 116, 127, 98, 105, 176, 187, 166, 173, 156, 151, 138, 129, 232, 227, 254, 245, 196, 207, 210, 217, 123, 112, 109, 102, 87, 92, 65, 74, 35, 40, 53, 62, 15, 4, 25, 18, 203, 192, 221, 214, 231, 236, 241, 250, 147, 152, 133, 142, 191, 180, 169, 162, 246, 253, 224, 235, 218, 209, 204, 199, 174, 165, 184, 179, 130, 137, 148, 159, 70, 77, 80, 91, 106, 97, 124, 119, 30, 21, 8, 3, 50, 57, 36, 47, 141, 134, 155, 144, 161, 170, 183, 188, 213, 222, 195, 200, 249, 242, 239, 228, 61, 54, 43, 32, 17, 26, 7, 12, 101, 110, 115, 120, 73, 66, 95, 84, 247, 252, 225, 234, 219, 208, 205, 198, 175, 164, 185, 178, 131, 136, 149, 158, 71, 76, 81, 90, 107, 96, 125, 118, 31, 20, 9, 2, 51, 56, 37, 46, 140, 135, 154, 145, 160, 171, 182, 189, 212, 223, 194, 201, 248, 243, 238, 229, 60, 55, 42, 33, 16, 27, 6, 13, 100, 111, 114, 121, 72, 67, 94, 85, 1, 10, 23, 28, 45, 38, 59, 48, 89, 82, 79, 68, 117, 126, 99, 104, 177, 186, 167, 172, 157, 150, 139, 128, 233, 226, 255, 244, 197, 206, 211, 216, 122, 113, 108, 103, 86, 93, 64, 75, 34, 41, 52, 63, 14, 5, 24, 19, 202, 193, 220, 215, 230, 237, 240, 251, 146, 153, 132, 143, 190, 181, 168, 163, 0, 13, 26, 23, 52, 57, 46, 35, 104, 101, 114, 127, 92, 81, 70, 75, 208, 221, 202, 199, 228, 233, 254, 243, 184, 181, 162, 175, 140, 129, 150, 155, 187, 182, 161, 172, 143, 130, 149, 152, 211, 222, 201, 196, 231, 234, 253, 240, 107, 102, 113, 124, 95, 82, 69, 72, 3, 14, 25, 20, 55, 58, 45, 32, 109, 96, 119, 122, 89, 84, 67, 78, 5, 8, 31, 18, 49, 60, 43, 38, 189, 176, 167, 170, 137, 132, 147, 158, 213, 216, 207, 194, 225, 236, 251, 246, 214, 219, 204, 193, 226, 239, 248, 245, 190, 179, 164, 169, 138, 135, 144, 157, 6, 11, 28, 17, 50, 63, 40, 37, 110, 99, 116, 121, 90, 87, 64, 77, 218, 215, 192, 205, 238, 227, 244, 249, 178, 191, 168, 165, 134, 139, 156, 145, 10, 7, 16, 29, 62, 51, 36, 41, 98, 111, 120, 117, 86, 91, 76, 65, 97, 108, 123, 118, 85, 88, 79, 66, 9, 4, 19, 30, 61, 48, 39, 42, 177, 188, 171, 166, 133, 136, 159, 146, 217, 212, 195, 206, 237, 224, 247, 250, 183, 186, 173, 160, 131, 142, 153, 148, 223, 210, 197, 200, 235, 230, 241, 252, 103, 106, 125, 112, 83, 94, 73, 68, 15, 2, 21, 24, 59, 54, 33, 44, 12, 1, 22, 27, 56, 53, 34, 47, 100, 105, 126, 115, 80, 93, 74, 71, 220, 209, 198, 203, 232, 229, 242, 255, 180, 185, 174, 163, 128, 141, 154, 151, 0, 14, 28, 18, 56, 54, 36, 42, 112, 126, 108, 98, 72, 70, 84, 90, 224, 238, 252, 242, 216, 214, 196, 202, 144, 158, 140, 130, 168, 166, 180, 186, 219, 213, 199, 201, 227, 237, 255, 241, 171, 165, 183, 185, 147, 157, 143, 129, 59, 53, 39, 41, 3, 13, 31, 17, 75, 69, 87, 89, 115, 125, 111, 97, 173, 163, 177, 191, 149, 155, 137, 135, 221, 211, 193, 207, 229, 235, 249, 247, 77, 67, 81, 95, 117, 123, 105, 103, 61, 51, 33, 47, 5, 11, 25, 23, 118, 120, 106, 100, 78, 64, 82, 92, 6, 8, 26, 20, 62, 48, 34, 44, 150, 152, 138, 132, 174, 160, 178, 188, 230, 232, 250, 244, 222, 208, 194, 204, 65, 79, 93, 83, 121, 119, 101, 107, 49, 63, 45, 35, 9, 7, 21, 27, 161, 175, 189, 179, 153, 151, 133, 139, 209, 223, 205, 195, 233, 231, 245, 251, 154, 148, 134, 136, 162, 172, 190, 176, 234, 228, 246, 248, 210, 220, 206, 192, 122, 116, 102, 104, 66, 76, 94, 80, 10, 4, 22, 24, 50, 60, 46, 32, 236, 226, 240, 254, 212, 218, 200, 198, 156, 146, 128, 142, 164, 170, 184, 182, 12, 2, 16, 30, 52, 58, 40, 38, 124, 114, 96, 110, 68, 74, 88, 86, 55, 57, 43, 37, 15, 1, 19, 29, 71, 73, 91, 85, 127, 113, 99, 109, 215, 217, 203, 197, 239, 225, 243, 253, 167, 169, 187, 181, 159, 145, 131, 141 ], kt = 2048, Ut = 16, xt = function(t) {
        m.call(this, t);
    }.prototype;
    xt.reset = S, xt.process = A, xt.finish = E;
    var Lt = function(t) {
        m.call(this, t);
    }.prototype;
    Lt.reset = S, Lt.process = _, Lt.finish = k;
    var Ht = m.prototype;
    Ht.reset = S, Ht.encrypt = function(t) {
        var e, r = A.call(this, t).result, i = E.call(this).result;
        return (e = new Uint8Array(r.length + i.length)).set(r), i.length > 0 && e.set(i, r.length), 
        this.result = e, this;
    }, Ht.decrypt = function(t) {
        var e, r = _.call(this, t).result, i = k.call(this).result;
        return (e = new Uint8Array(r.length + i.length)).set(r), i.length > 0 && e.set(i, r.length), 
        this.result = e, this;
    };
    var It = 68719476704, qt = U.prototype;
    qt.reset = L, qt.encrypt = function(t) {
        var e = H.call(this, t).result, r = I.call(this).result, i = new Uint8Array(e.length + r.length);
        return e.length && i.set(e), r.length && i.set(r, e.length), this.result = i, this;
    }, qt.decrypt = function(t) {
        var e = q.call(this, t).result, r = M.call(this).result, i = new Uint8Array(e.length + r.length);
        return e.length && i.set(e), r.length && i.set(r, e.length), this.result = i, this;
    };
    var Mt = function(t) {
        U.call(this, t);
    }.prototype;
    Mt.reset = L, Mt.process = H, Mt.finish = I;
    var Ct = function(t) {
        U.call(this, t);
    }.prototype;
    Ct.reset = L, Ct.process = q, Ct.finish = M;
    var p, Zt = new Uint8Array(1048576), zt = new m({
        heap: Zt,
        asm: p = v(e, null, Zt.buffer)
    });
    t.AES_CBC = {
        encrypt: function(t, e, r, i) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            return zt.reset({
                key: e,
                padding: r,
                iv: i
            }).encrypt(t).result;
        },
        decrypt: function(t, e, r, i) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            return zt.reset({
                key: e,
                padding: r,
                iv: i
            }).decrypt(t).result;
        }
    };
    var Ot = new U({
        heap: Zt,
        asm: p
    });
    t.AES_GCM = {
        encrypt: function(t, e, r, i, n) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            if (void 0 === r) throw new SyntaxError("iv required");
            return Ot.reset({
                key: e,
                iv: r,
                adata: i,
                tagSize: n
            }).encrypt(t).result;
        },
        decrypt: function(t, e, r, i, n) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            if (void 0 === r) throw new SyntaxError("iv required");
            return Ot.reset({
                key: e,
                iv: r,
                adata: i,
                tagSize: n
            }).decrypt(t).result;
        }
    };
    var Tt = 64, Rt = 20;
    C.BLOCK_SIZE = Tt, C.HASH_SIZE = Rt;
    var Bt = C.prototype;
    Bt.reset = function() {
        return this.result = null, this.pos = 0, this.len = 0, this.asm.reset(), this;
    }, Bt.process = function(t) {
        if (null !== this.result) throw new r("state must be reset before processing new data");
        if (l(t) && (t = s(t)), c(t) && (t = new Uint8Array(t)), !w(t)) throw new TypeError("data isn't of expected type");
        for (var e = this.asm, i = this.heap, n = this.pos, a = this.len, h = 0, o = t.length, f = 0; o > 0; ) f = o > (f = i.length - n - a) ? f : o, 
        i.set(new Uint8Array(t.buffer || t, h, f), this.pos + this.len), i.set(t.subarray(h, h + f), n + a), 
        a += f, h += f, o -= f, n += f = e.process(n, a), (a -= f) || (n = 0);
        return this.pos = n, this.len = a, this;
    }, Bt.finish = function() {
        if (null !== this.result) throw new r("state must be reset before processing new data");
        return this.asm.finish(this.pos, this.len, 0), this.result = new Uint8Array(Rt), 
        this.result.set(this.heap.subarray(0, Rt)), this.pos = 0, this.len = 0, this;
    };
    var Pt = null;
    t.SHA1 = {
        bytes: z,
        hex: function(t) {
            return h(z(t));
        },
        base64: function(t) {
            return o(z(t));
        }
    };
    var Kt = 64, Nt = 32;
    O.BLOCK_SIZE = Kt, O.HASH_SIZE = Nt;
    var jt = O.prototype;
    jt.reset = function() {
        return this.result = null, this.pos = 0, this.len = 0, this.asm.reset(), this;
    }, jt.process = function(t) {
        if (null !== this.result) throw new r("state must be reset before processing new data");
        if (l(t) && (t = s(t)), c(t) && (t = new Uint8Array(t)), !w(t)) throw new TypeError("data isn't of expected type");
        for (var e = this.asm, i = this.heap, n = this.pos, a = this.len, h = 0, o = t.length, f = 0; o > 0; ) f = o > (f = i.length - n - a) ? f : o, 
        i.set(new Uint8Array(t.buffer || t, h, f), this.pos + this.len), i.set(t.subarray(h, h + f), n + a), 
        a += f, h += f, o -= f, n += f = e.process(n, a), (a -= f) || (n = 0);
        return this.pos = n, this.len = a, this;
    }, jt.finish = function() {
        if (null !== this.result) throw new r("state must be reset before processing new data");
        return this.asm.finish(this.pos, this.len, 0), this.result = new Uint8Array(Nt), 
        this.result.set(this.heap.subarray(0, Nt)), this.pos = 0, this.len = 0, this;
    };
    var Dt = null;
    t.SHA256 = {
        bytes: R,
        hex: function(t) {
            return h(R(t));
        },
        base64: function(t) {
            return o(R(t));
        }
    };
    var Gt = B.prototype;
    Gt.reset = function(t) {
        var e = (t = t || {}).password;
        if (null === this.key && !l(e) && !e) throw new r("no key is associated with the instance");
        this.result = null, this.hash.reset(), (e || l(e)) && (this.key = P(this.hash, e));
        for (var i = new Uint8Array(this.key), n = 0; n < i.length; ++n) i[n] ^= 54;
        this.hash.process(i);
        var s = t.verify;
        return void 0 !== s ? K.call(this, s) : this.verify = null, this;
    }, Gt.process = N, Gt.finish = function() {
        if (null === this.key) throw new r("no key is associated with the instance");
        if (null !== this.result) throw new r("state must be reset before processing new data");
        for (var t = this.hash.finish().result, e = new Uint8Array(this.key), i = 0; i < e.length; ++i) e[i] ^= 92;
        var n = this.verify, s = this.hash.reset().process(e).process(t).finish().result;
        if (n) if (n.length === s.length) {
            var a = 0;
            for (i = 0; i < n.length; i++) a |= n[i] ^ s[i];
            this.result = !a;
        } else this.result = !1; else this.result = s;
        return this;
    }, j.BLOCK_SIZE = C.BLOCK_SIZE, j.HMAC_SIZE = C.HASH_SIZE;
    var Ft = j.prototype;
    Ft.reset = function(t) {
        t = t || {}, this.result = null, this.hash.reset();
        var e = t.password;
        if (void 0 !== e) {
            l(e) && (e = s(e));
            var r = this.key = P(this.hash, e);
            this.hash.reset().asm.hmac_init(r[0] << 24 | r[1] << 16 | r[2] << 8 | r[3], r[4] << 24 | r[5] << 16 | r[6] << 8 | r[7], r[8] << 24 | r[9] << 16 | r[10] << 8 | r[11], r[12] << 24 | r[13] << 16 | r[14] << 8 | r[15], r[16] << 24 | r[17] << 16 | r[18] << 8 | r[19], r[20] << 24 | r[21] << 16 | r[22] << 8 | r[23], r[24] << 24 | r[25] << 16 | r[26] << 8 | r[27], r[28] << 24 | r[29] << 16 | r[30] << 8 | r[31], r[32] << 24 | r[33] << 16 | r[34] << 8 | r[35], r[36] << 24 | r[37] << 16 | r[38] << 8 | r[39], r[40] << 24 | r[41] << 16 | r[42] << 8 | r[43], r[44] << 24 | r[45] << 16 | r[46] << 8 | r[47], r[48] << 24 | r[49] << 16 | r[50] << 8 | r[51], r[52] << 24 | r[53] << 16 | r[54] << 8 | r[55], r[56] << 24 | r[57] << 16 | r[58] << 8 | r[59], r[60] << 24 | r[61] << 16 | r[62] << 8 | r[63]);
        } else this.hash.asm.hmac_reset();
        var i = t.verify;
        return void 0 !== i ? K.call(this, i) : this.verify = null, this;
    }, Ft.process = N, Ft.finish = function() {
        if (null === this.key) throw new r("no key is associated with the instance");
        if (null !== this.result) throw new r("state must be reset before processing new data");
        var t = this.hash, e = this.hash.asm, i = this.hash.heap;
        e.hmac_finish(t.pos, t.len, 0);
        var n = this.verify, s = new Uint8Array(Rt);
        if (s.set(i.subarray(0, Rt)), n) if (n.length === s.length) {
            for (var a = 0, h = 0; h < n.length; h++) a |= n[h] ^ s[h];
            this.result = !a;
        } else this.result = !1; else this.result = s;
        return this;
    };
    var Wt = null;
    G.BLOCK_SIZE = O.BLOCK_SIZE, G.HMAC_SIZE = O.HASH_SIZE;
    var Vt = G.prototype;
    Vt.reset = function(t) {
        t = t || {}, this.result = null, this.hash.reset();
        var e = t.password;
        if (void 0 !== e) {
            l(e) && (e = s(e));
            var r = this.key = P(this.hash, e);
            this.hash.reset().asm.hmac_init(r[0] << 24 | r[1] << 16 | r[2] << 8 | r[3], r[4] << 24 | r[5] << 16 | r[6] << 8 | r[7], r[8] << 24 | r[9] << 16 | r[10] << 8 | r[11], r[12] << 24 | r[13] << 16 | r[14] << 8 | r[15], r[16] << 24 | r[17] << 16 | r[18] << 8 | r[19], r[20] << 24 | r[21] << 16 | r[22] << 8 | r[23], r[24] << 24 | r[25] << 16 | r[26] << 8 | r[27], r[28] << 24 | r[29] << 16 | r[30] << 8 | r[31], r[32] << 24 | r[33] << 16 | r[34] << 8 | r[35], r[36] << 24 | r[37] << 16 | r[38] << 8 | r[39], r[40] << 24 | r[41] << 16 | r[42] << 8 | r[43], r[44] << 24 | r[45] << 16 | r[46] << 8 | r[47], r[48] << 24 | r[49] << 16 | r[50] << 8 | r[51], r[52] << 24 | r[53] << 16 | r[54] << 8 | r[55], r[56] << 24 | r[57] << 16 | r[58] << 8 | r[59], r[60] << 24 | r[61] << 16 | r[62] << 8 | r[63]);
        } else this.hash.asm.hmac_reset();
        var i = t.verify;
        return void 0 !== i ? K.call(this, i) : this.verify = null, this;
    }, Vt.process = N, Vt.finish = function() {
        if (null === this.key) throw new r("no key is associated with the instance");
        if (null !== this.result) throw new r("state must be reset before processing new data");
        var t = this.hash, e = this.hash.asm, i = this.hash.heap;
        e.hmac_finish(t.pos, t.len, 0);
        var n = this.verify, s = new Uint8Array(Nt);
        if (s.set(i.subarray(0, Nt)), n) if (n.length === s.length) {
            for (var a = 0, h = 0; h < n.length; h++) a |= n[h] ^ s[h];
            this.result = !a;
        } else this.result = !1; else this.result = s;
        return this;
    };
    var Jt = null;
    t.HMAC = t.HMAC_SHA1 = {
        bytes: W,
        hex: function(t, e) {
            return h(W(t, e));
        },
        base64: function(t, e) {
            return o(W(t, e));
        }
    }, t.HMAC_SHA256 = {
        bytes: V,
        hex: function(t, e) {
            return h(V(t, e));
        },
        base64: function(t, e) {
            return o(V(t, e));
        }
    };
    var Qt = J.prototype;
    Qt.reset = Q, Qt.generate = function(t, e, n) {
        if (null !== this.result) throw new r("state must be reset before processing new data");
        if (!t && !l(t)) throw new i("bad 'salt' value");
        e = e || this.count, n = n || this.length, this.result = new Uint8Array(n);
        for (var s = Math.ceil(n / this.hmac.HMAC_SIZE), a = 1; s >= a; ++a) {
            var h = (a - 1) * this.hmac.HMAC_SIZE, o = (s > a ? 0 : n % this.hmac.HMAC_SIZE) || this.hmac.HMAC_SIZE, f = new Uint8Array(this.hmac.reset().process(t).process(new Uint8Array([ a >>> 24 & 255, a >>> 16 & 255, a >>> 8 & 255, 255 & a ])).finish().result);
            this.result.set(f.subarray(0, o), h);
            for (var u = 1; e > u; ++u) {
                f = new Uint8Array(this.hmac.reset().process(f).finish().result);
                for (var c = 0; o > c; ++c) this.result[h + c] ^= f[c];
            }
        }
        return this;
    };
    var Xt = X.prototype;
    Xt.reset = Q, Xt.generate = function(t, e, n) {
        if (null !== this.result) throw new r("state must be reset before processing new data");
        if (!t && !l(t)) throw new i("bad 'salt' value");
        e = e || this.count, n = n || this.length, this.result = new Uint8Array(n);
        for (var s = Math.ceil(n / this.hmac.HMAC_SIZE), a = 1; s >= a; ++a) {
            var h = (a - 1) * this.hmac.HMAC_SIZE, o = (s > a ? 0 : n % this.hmac.HMAC_SIZE) || this.hmac.HMAC_SIZE;
            this.hmac.reset().process(t), this.hmac.hash.asm.pbkdf2_generate_block(this.hmac.hash.pos, this.hmac.hash.len, a, e, 0), 
            this.result.set(this.hmac.hash.heap.subarray(0, o), h);
        }
        return this;
    };
    var Yt = null, $t = Y.prototype;
    $t.reset = Q, $t.generate = function(t, e, n) {
        if (null !== this.result) throw new r("state must be reset before processing new data");
        if (!t && !l(t)) throw new i("bad 'salt' value");
        e = e || this.count, n = n || this.length, this.result = new Uint8Array(n);
        for (var s = Math.ceil(n / this.hmac.HMAC_SIZE), a = 1; s >= a; ++a) {
            var h = (a - 1) * this.hmac.HMAC_SIZE, o = (s > a ? 0 : n % this.hmac.HMAC_SIZE) || this.hmac.HMAC_SIZE;
            this.hmac.reset().process(t), this.hmac.hash.asm.pbkdf2_generate_block(this.hmac.hash.pos, this.hmac.hash.len, a, e, 0), 
            this.result.set(this.hmac.hash.heap.subarray(0, o), h);
        }
        return this;
    };
    var te = null;
    t.PBKDF2 = t.PBKDF2_HMAC_SHA1 = {
        bytes: tt,
        hex: function(t, e, r, i) {
            return h(tt(t, e, r, i));
        },
        base64: function(t, e, r, i) {
            return o(tt(t, e, r, i));
        }
    }, t.PBKDF2_HMAC_SHA256 = {
        bytes: et,
        hex: function(t, e, r, i) {
            return h(et(t, e, r, i));
        },
        base64: function(t, e, r, i) {
            return o(et(t, e, r, i));
        }
    };
    var ee, re = function() {
        function t() {
            function t() {
                s ^= f << 11, f = f + u | 0, f ^= u >>> 2, u = u + (l = l + s | 0) | 0, u ^= l << 8, 
                l = l + (c = c + f | 0) | 0, l ^= c >>> 16, c = c + (w = w + u | 0) | 0, c ^= w << 10, 
                w = w + (y = y + l | 0) | 0, w ^= y >>> 4, y = y + (p = p + c | 0) | 0, y ^= p << 8, 
                p = p + (s = s + w | 0) | 0, u = u + (p ^= s >>> 9) | 0, s = s + (f = f + y | 0) | 0;
            }
            var s, f, u, l, c, w, y, p;
            n = a = h = 0, s = f = u = l = c = w = y = p = 2654435769;
            for (var v = 0; 4 > v; v++) t();
            for (v = 0; 256 > v; v += 8) s = s + i[0 | v] | 0, f = f + i[1 | v] | 0, u = u + i[2 | v] | 0, 
            l = l + i[3 | v] | 0, c = c + i[4 | v] | 0, w = w + i[5 | v] | 0, y = y + i[6 | v] | 0, 
            p = p + i[7 | v] | 0, t(), r.set([ s, f, u, l, c, w, y, p ], v);
            for (v = 0; 256 > v; v += 8) s = s + r[0 | v] | 0, f = f + r[1 | v] | 0, u = u + r[2 | v] | 0, 
            l = l + r[3 | v] | 0, c = c + r[4 | v] | 0, w = w + r[5 | v] | 0, y = y + r[6 | v] | 0, 
            p = p + r[7 | v] | 0, t(), r.set([ s, f, u, l, c, w, y, p ], v);
            e(1), o = 256;
        }
        function e(t) {
            t = t || 1;
            for (var e, s, o; t--; ) for (a = a + (h = h + 1 | 0) | 0, e = 0; 256 > e; e += 4) n ^= n << 13, 
            n = r[e + 128 & 255] + n | 0, s = r[0 | e], r[0 | e] = o = r[s >>> 2 & 255] + (n + a | 0) | 0, 
            i[0 | e] = a = r[o >>> 10 & 255] + s | 0, n ^= n >>> 6, n = r[e + 129 & 255] + n | 0, 
            s = r[1 | e], r[1 | e] = o = r[s >>> 2 & 255] + (n + a | 0) | 0, i[1 | e] = a = r[o >>> 10 & 255] + s | 0, 
            n ^= n << 2, n = r[e + 130 & 255] + n | 0, s = r[2 | e], r[2 | e] = o = r[s >>> 2 & 255] + (n + a | 0) | 0, 
            i[2 | e] = a = r[o >>> 10 & 255] + s | 0, n ^= n >>> 16, n = r[e + 131 & 255] + n | 0, 
            s = r[3 | e], r[3 | e] = o = r[s >>> 2 & 255] + (n + a | 0) | 0, i[3 | e] = a = r[o >>> 10 & 255] + s | 0;
        }
        var r = new Uint32Array(256), i = new Uint32Array(256), n = 0, a = 0, h = 0, o = 0;
        return {
            seed: function(e) {
                var r, n, a, h, o;
                if (y(e)) e = new Uint8Array(e.buffer); else if (u(e)) (h = new Et(1))[0] = e, e = new Uint8Array(h.buffer); else if (l(e)) e = s(e); else {
                    if (!c(e)) throw new TypeError("bad seed type");
                    e = new Uint8Array(e);
                }
                for (o = e.length, n = 0; o > n; n += 1024) {
                    for (a = n, r = 0; 1024 > r && o > a; a = n | ++r) i[r >> 2] ^= e[a] << ((3 & r) << 3);
                    t();
                }
            },
            prng: e,
            rand: function() {
                return o-- || (e(1), o = 255), i[o];
            }
        };
    }(), ie = e.console, ne = e.Date.now, se = e.Math.random, ae = e.performance, he = e.crypto || e.msCrypto;
    void 0 !== he && (ee = he.getRandomValues);
    var oe, fe, ue = re.rand, le = re.seed, ce = 0, we = !1, ye = !1, pe = 0, ve = 256, de = !1, ge = !1, be = {};
    if (void 0 !== ae) oe = function() {
        return 1e3 * ae.now() | 0;
    }; else {
        var me = 1e3 * ne() | 0;
        oe = function() {
            return 1e3 * ne() - me | 0;
        };
    }
    t.random = st, t.random.seed = it, Object.defineProperty(st, "allowWeak", {
        get: function() {
            return de;
        },
        set: function(t) {
            de = t;
        }
    }), Object.defineProperty(st, "skipSystemRNGWarning", {
        get: function() {
            return ge;
        },
        set: function(t) {
            ge = t;
        }
    }), t.getRandomValues = nt, t.getRandomValues.seed = it, Object.defineProperty(nt, "allowWeak", {
        get: function() {
            return de;
        },
        set: function(t) {
            de = t;
        }
    }), Object.defineProperty(nt, "skipSystemRNGWarning", {
        get: function() {
            return ge;
        },
        set: function(t) {
            ge = t;
        }
    }), e.Math.random = st, void 0 === e.crypto && (e.crypto = {}), e.crypto.getRandomValues = nt, 
    fe = void 0 === e.Math.imul ? function(t, r, i) {
        e.Math.imul = at;
        var n = ht(t, r, i);
        return delete e.Math.imul, n;
    } : ht;
    var Se = new Uint32Array(1048576), ht = fe(e, null, Se.buffer), Ae = new Uint32Array(0), Ee = ft.prototype = new Number();
    Ee.toString = function(t) {
        t = t || 16;
        var e = this.limbs, r = this.bitLength, n = "";
        if (16 !== t) throw new i("bad radix");
        for (var s = (r + 31 >> 5) - 1; s >= 0; s--) {
            var a = e[s].toString(16);
            n += "00000000".substr(a.length), n += a;
        }
        return (n = n.replace(/^0+/, "")).length || (n = "0"), this.sign < 0 && (n = "-" + n), 
        n;
    }, Ee.toBytes = function() {
        var t = this.bitLength, e = this.limbs;
        if (0 === t) return new Uint8Array(0);
        for (var r = t + 7 >> 3, i = new Uint8Array(r), n = 0; r > n; n++) {
            var s = r - n - 1;
            i[n] = e[s >> 2] >> ((3 & s) << 3);
        }
        return i;
    }, Ee.valueOf = function() {
        var t = this.limbs, e = this.bitLength, r = this.sign;
        if (!r) return 0;
        if (32 >= e) return r * (t[0] >>> 0);
        if (52 >= e) return r * (4294967296 * (t[1] >>> 0) + (t[0] >>> 0));
        var i, n, s = 0;
        for (i = t.length - 1; i >= 0; i--) if (0 !== (n = t[i])) {
            for (;0 == (n << s & 2147483648); ) s++;
            break;
        }
        return 0 === i ? r * (t[0] >>> 0) : r * (1048576 * ((t[i] << s | (s ? t[i - 1] >>> 32 - s : 0)) >>> 0) + ((t[i - 1] << s | (s && i > 1 ? t[i - 2] >>> 32 - s : 0)) >>> 12)) * Math.pow(2, 32 * i - s - 52);
    }, Ee.clamp = function(t) {
        var e = this.limbs;
        if (t >= this.bitLength) return this;
        var r = new ft(), i = t + 31 >> 5, n = t % 32;
        return r.limbs = new Uint32Array(e.subarray(0, i)), r.bitLength = t, r.sign = this.sign, 
        n && (r.limbs[i - 1] &= -1 >>> 32 - n), r;
    }, Ee.slice = function(t, e) {
        if (!u(t)) throw new TypeError("TODO");
        if (void 0 !== e && !u(e)) throw new TypeError("TODO");
        var r = this.limbs, i = this.bitLength;
        if (0 > t) throw new RangeError("TODO");
        if (t >= i) return _e;
        (void 0 === e || e > i - t) && (e = i - t);
        var n, s = new ft(), a = t >> 5, h = t + e + 31 >> 5, o = e + 31 >> 5, f = t % 32, l = e % 32;
        if (n = new Uint32Array(o), f) {
            for (var c = 0; h - a - 1 > c; c++) n[c] = r[a + c] >>> f | r[a + c + 1] << 32 - f;
            n[c] = r[a + c] >>> f;
        } else n.set(r.subarray(a, h));
        return l && (n[o - 1] &= -1 >>> 32 - l), s.limbs = n, s.bitLength = e, s.sign = this.sign, 
        s;
    }, Ee.negate = function() {
        var t = new ft();
        return t.limbs = this.limbs, t.bitLength = this.bitLength, t.sign = -1 * this.sign, 
        t;
    }, Ee.compare = function(t) {
        ot(t) || (t = new ft(t));
        var e = this.limbs, r = e.length, i = t.limbs, n = i.length;
        return this.sign < t.sign ? -1 : this.sign > t.sign ? 1 : (Se.set(e, 0), Se.set(i, r), 
        ht.cmp(0, r << 2, r << 2, n << 2) * this.sign);
    }, Ee.add = function(t) {
        if (ot(t) || (t = new ft(t)), !this.sign) return t;
        if (!t.sign) return this;
        var e, r, i, n, s = this.bitLength, a = this.limbs, h = a.length, o = this.sign, f = t.bitLength, u = t.limbs, l = u.length, c = t.sign, w = new ft();
        r = (e = (s > f ? s : f) + (o * c > 0 ? 1 : 0)) + 31 >> 5, ht.sreset();
        var y = ht.salloc(h << 2), p = ht.salloc(l << 2), v = ht.salloc(r << 2);
        return ht.z(v - y + (r << 2), 0, y), Se.set(a, y >> 2), Se.set(u, p >> 2), o * c > 0 ? (ht.add(y, h << 2, p, l << 2, v, r << 2), 
        i = o) : i = o > c ? (n = ht.sub(y, h << 2, p, l << 2, v, r << 2)) ? c : o : (n = ht.sub(p, l << 2, y, h << 2, v, r << 2)) ? o : c, 
        n && ht.neg(v, r << 2, v, r << 2), 0 === ht.tst(v, r << 2) ? _e : (w.limbs = new Uint32Array(Se.subarray(v >> 2, (v >> 2) + r)), 
        w.bitLength = e, w.sign = i, w);
    }, Ee.subtract = function(t) {
        return ot(t) || (t = new ft(t)), this.add(t.negate());
    }, Ee.multiply = function(t) {
        if (ot(t) || (t = new ft(t)), !this.sign || !t.sign) return _e;
        var e, r, i = this.bitLength, n = this.limbs, s = n.length, a = t.bitLength, h = t.limbs, o = h.length, f = new ft();
        r = (e = i + a) + 31 >> 5, ht.sreset();
        var u = ht.salloc(s << 2), l = ht.salloc(o << 2), c = ht.salloc(r << 2);
        return ht.z(c - u + (r << 2), 0, u), Se.set(n, u >> 2), Se.set(h, l >> 2), ht.mul(u, s << 2, l, o << 2, c, r << 2), 
        f.limbs = new Uint32Array(Se.subarray(c >> 2, (c >> 2) + r)), f.sign = this.sign * t.sign, 
        f.bitLength = e, f;
    }, Ee.square = function() {
        if (!this.sign) return _e;
        var t, e, r = this.bitLength, i = this.limbs, n = i.length, s = new ft();
        e = 31 + (t = r << 1) >> 5, ht.sreset();
        var a = ht.salloc(n << 2), h = ht.salloc(e << 2);
        return ht.z(h - a + (e << 2), 0, a), Se.set(i, a >> 2), ht.sqr(a, n << 2, h), s.limbs = new Uint32Array(Se.subarray(h >> 2, (h >> 2) + e)), 
        s.bitLength = t, s.sign = 1, s;
    }, Ee.divide = function(t) {
        ot(t) || (t = new ft(t));
        var e, r, i = this.bitLength, n = this.limbs, s = n.length, a = t.bitLength, h = t.limbs, o = h.length, f = _e, u = _e;
        ht.sreset();
        var l = ht.salloc(s << 2), c = ht.salloc(o << 2), w = ht.salloc(o << 2), y = ht.salloc(s << 2);
        return ht.z(y - l + (s << 2), 0, l), Se.set(n, l >> 2), Se.set(h, c >> 2), ht.div(l, s << 2, c, o << 2, w, y), 
        (e = ht.tst(y, s << 2) >> 2) && ((f = new ft()).limbs = new Uint32Array(Se.subarray(y >> 2, (y >> 2) + e)), 
        f.bitLength = e << 5 > i ? i : e << 5, f.sign = this.sign * t.sign), (r = ht.tst(w, o << 2) >> 2) && ((u = new ft()).limbs = new Uint32Array(Se.subarray(w >> 2, (w >> 2) + r)), 
        u.bitLength = r << 5 > a ? a : r << 5, u.sign = this.sign), {
            quotient: f,
            remainder: u
        };
    };
    var _e = new ft(0), ke = new ft(1);
    Object.freeze(_e), Object.freeze(ke);
    var Ue = lt.prototype = new ft();
    Ue.reduce = function(t) {
        return ot(t) || (t = new ft(t)), t.bitLength <= 32 && this.bitLength <= 32 ? new ft(t.valueOf() % this.valueOf()) : t.compare(this) < 0 ? t : t.divide(this).remainder;
    }, Ue.inverse = function(t) {
        var e = ut(this, t = this.reduce(t));
        return 1 !== e.gcd.valueOf() ? null : ((e = e.y).sign < 0 && (e = e.add(this).clamp(this.bitLength)), 
        e);
    }, Ue.power = function(t, e) {
        ot(t) || (t = new ft(t)), ot(e) || (e = new ft(e));
        for (var r = 0, i = 0; i < e.limbs.length; i++) for (var n = e.limbs[i]; n; ) 1 & n && r++, 
        n >>>= 1;
        var s = 8;
        e.bitLength <= 4536 && (s = 7), e.bitLength <= 1736 && (s = 6), e.bitLength <= 630 && (s = 5), 
        e.bitLength <= 210 && (s = 4), e.bitLength <= 60 && (s = 3), e.bitLength <= 12 && (s = 2), 
        1 << s - 1 >= r && (s = 1);
        var a = ct((t = ct(this.reduce(t).multiply(this.comodulusRemainderSquare), this)).square(), this), h = new Array(1 << s - 1);
        for (h[0] = t, h[1] = ct(t.multiply(a), this), i = 2; 1 << s - 1 > i; i++) h[i] = ct(h[i - 1].multiply(a), this);
        var o = this.comodulusRemainder, f = o;
        for (i = e.limbs.length - 1; i >= 0; i--) {
            n = e.limbs[i];
            for (var u = 32; u > 0; ) if (2147483648 & n) {
                for (var l = n >>> 32 - s, c = s; 0 == (1 & l); ) l >>>= 1, c--;
                for (var w = h[l >>> 1]; l; ) l >>>= 1, f !== o && (f = ct(f.square(), this));
                f = f !== o ? ct(f.multiply(w), this) : w, n <<= c, u -= c;
            } else f !== o && (f = ct(f.square(), this)), n <<= 1, u--;
        }
        return ct(f, this);
    };
    var xe = [ 2, 3 ];
    Ee.isProbablePrime = function(t) {
        t = t || 80;
        var e = this.limbs, r = 0;
        if (0 == (1 & e[0])) return !1;
        if (1 >= t) return !0;
        var i = 0, n = 0, s = 0;
        for (r = 0; r < e.length; r++) {
            for (var a = e[r]; a; ) i += 3 & a, a >>>= 2;
            for (var h = e[r]; h; ) n += 3 & h, n -= 3 & (h >>>= 2), h >>>= 2;
            for (var o = e[r]; o; ) s += 15 & o, s -= 15 & (o >>>= 4), o >>>= 4;
        }
        return !!(i % 3 && n % 5 && s % 17) && (2 >= t || wt.call(this, t >>> 1));
    }, ft.randomProbablePrime = yt, ft.ZERO = _e, ft.ONE = ke, ft.extGCD = ut, t.BigNumber = ft, 
    t.Modulus = lt;
    var Le = pt.prototype;
    Le.reset = vt, Le.encrypt = dt, Le.decrypt = gt, pt.generateKey = bt;
    var He = mt.prototype;
    He.reset = function(t) {
        var e = (t = t || {}).label;
        if (void 0 !== e) {
            if (c(e) || w(e)) e = new Uint8Array(e); else {
                if (!l(e)) throw new TypeError("unexpected label type");
                e = s(e);
            }
            this.label = e.length > 0 ? e : null;
        } else this.label = null;
        vt.call(this, t);
    }, He.encrypt = function(t) {
        if (!this.key) throw new r("no key is associated with the instance");
        var e = Math.ceil(this.key[0].bitLength / 8), n = this.hash.HASH_SIZE, a = t.byteLength || t.length || 0, h = e - a - 2 * n - 2;
        if (a > e - 2 * this.hash.HASH_SIZE - 2) throw new i("data too large");
        var o = new Uint8Array(e), f = o.subarray(1, n + 1), u = o.subarray(n + 1);
        if (w(t)) u.set(t, n + h + 1); else if (c(t)) u.set(new Uint8Array(t), n + h + 1); else {
            if (!l(t)) throw new TypeError("unexpected data type");
            u.set(s(t), n + h + 1);
        }
        u.set(this.hash.reset().process(this.label || "").finish().result, 0), u[n + h] = 1, 
        nt(f);
        for (var y = St.call(this, f, u.length), p = 0; p < u.length; p++) u[p] ^= y[p];
        var v = St.call(this, u, f.length);
        for (p = 0; p < f.length; p++) f[p] ^= v[p];
        return dt.call(this, o), this;
    }, He.decrypt = function(t) {
        if (!this.key) throw new r("no key is associated with the instance");
        var e = Math.ceil(this.key[0].bitLength / 8), s = this.hash.HASH_SIZE;
        if ((t.byteLength || t.length || 0) !== e) throw new i("bad data");
        gt.call(this, t);
        var a = this.result[0], h = this.result.subarray(1, s + 1), o = this.result.subarray(s + 1);
        if (0 !== a) throw new n("decryption failed");
        for (var f = St.call(this, o, h.length), u = 0; u < h.length; u++) h[u] ^= f[u];
        var l = St.call(this, h, o.length);
        for (u = 0; u < o.length; u++) o[u] ^= l[u];
        var c = this.hash.reset().process(this.label || "").finish().result;
        for (u = 0; s > u; u++) if (c[u] !== o[u]) throw new n("decryption failed");
        for (var w = s; w < o.length; w++) {
            var y = o[w];
            if (1 === y) break;
            if (0 !== y) throw new n("decryption failed");
        }
        if (w === o.length) throw new n("decryption failed");
        return this.result = o.subarray(w + 1), this;
    };
    var Ie = At.prototype;
    Ie.reset = function(t) {
        t = t || {}, vt.call(this, t);
        var e = t.saltLength;
        if (void 0 !== e) {
            if (!u(e) || 0 > e) throw new TypeError("saltLength should be a non-negative number");
            if (null !== this.key && Math.ceil((this.key[0].bitLength - 1) / 8) < this.hash.HASH_SIZE + e + 2) throw new SyntaxError("saltLength is too large");
            this.saltLength = e;
        } else this.saltLength = 4;
    }, Ie.sign = function(t) {
        if (!this.key) throw new r("no key is associated with the instance");
        var e = this.key[0].bitLength, i = this.hash.HASH_SIZE, n = Math.ceil((e - 1) / 8), s = this.saltLength, a = n - s - i - 2, h = new Uint8Array(n), o = h.subarray(n - i - 1, n - 1), f = h.subarray(0, n - i - 1), u = f.subarray(a + 1), l = new Uint8Array(8 + i + s), c = l.subarray(8, 8 + i), w = l.subarray(8 + i);
        c.set(this.hash.reset().process(t).finish().result), s > 0 && nt(w), f[a] = 1, u.set(w), 
        o.set(this.hash.reset().process(l).finish().result);
        for (var y = St.call(this, o, f.length), p = 0; p < f.length; p++) f[p] ^= y[p];
        h[n - 1] = 188;
        var v = 8 * n - e + 1;
        return v % 8 && (h[0] &= 255 >>> v), gt.call(this, h), this;
    }, Ie.verify = function(t, e) {
        if (!this.key) throw new r("no key is associated with the instance");
        var i = this.key[0].bitLength, s = this.hash.HASH_SIZE, a = Math.ceil((i - 1) / 8), h = this.saltLength, o = a - h - s - 2;
        dt.call(this, t);
        var f = this.result;
        if (188 !== f[a - 1]) throw new n("bad signature");
        var u = f.subarray(a - s - 1, a - 1), l = f.subarray(0, a - s - 1), c = l.subarray(o + 1), w = 8 * a - i + 1;
        if (w % 8 && f[0] >>> 8 - w) throw new n("bad signature");
        for (var y = St.call(this, u, l.length), p = 0; p < l.length; p++) l[p] ^= y[p];
        for (w % 8 && (f[0] &= 255 >>> w), p = 0; o > p; p++) if (0 !== l[p]) throw new n("bad signature");
        if (1 !== l[o]) throw new n("bad signature");
        var v = new Uint8Array(8 + s + h), d = v.subarray(8, 8 + s), g = v.subarray(8 + s);
        d.set(this.hash.reset().process(e).finish().result), g.set(c);
        var b = this.hash.reset().process(v).finish().result;
        for (p = 0; s > p; p++) if (u[p] !== b[p]) throw new n("bad signature");
        return this;
    }, t.RSA = {
        generateKey: function(t, e) {
            if (void 0 === t) throw new SyntaxError("bitlen required");
            if (void 0 === e) throw new SyntaxError("e required");
            for (var r = bt(t, e), i = 0; i < r.length; i++) ot(r[i]) && (r[i] = r[i].toBytes());
            return r;
        }
    }, t.RSA_OAEP_SHA1 = {
        encrypt: function(t, e, r) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            return new mt({
                hash: Z(),
                key: e,
                label: r
            }).encrypt(t).result;
        },
        decrypt: function(t, e, r) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            return new mt({
                hash: Z(),
                key: e,
                label: r
            }).decrypt(t).result;
        }
    }, t.RSA_OAEP_SHA256 = {
        encrypt: function(t, e, r) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            return new mt({
                hash: T(),
                key: e,
                label: r
            }).encrypt(t).result;
        },
        decrypt: function(t, e, r) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            return new mt({
                hash: T(),
                key: e,
                label: r
            }).decrypt(t).result;
        }
    }, t.RSA_PSS_SHA1 = {
        sign: function(t, e, r) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            return new At({
                hash: Z(),
                key: e,
                saltLength: r
            }).sign(t).result;
        },
        verify: function(t, e, r, i) {
            if (void 0 === t) throw new SyntaxError("signature required");
            if (void 0 === e) throw new SyntaxError("data required");
            if (void 0 === r) throw new SyntaxError("key required");
            try {
                return new At({
                    hash: Z(),
                    key: r,
                    saltLength: i
                }).verify(t, e), !0;
            } catch (t) {
                if (!(t instanceof n)) throw t;
            }
            return !1;
        }
    }, t.RSA_PSS_SHA256 = {
        sign: function(t, e, r) {
            if (void 0 === t) throw new SyntaxError("data required");
            if (void 0 === e) throw new SyntaxError("key required");
            return new At({
                hash: T(),
                key: e,
                saltLength: r
            }).sign(t).result;
        },
        verify: function(t, e, r, i) {
            if (void 0 === t) throw new SyntaxError("signature required");
            if (void 0 === e) throw new SyntaxError("data required");
            if (void 0 === r) throw new SyntaxError("key required");
            try {
                return new At({
                    hash: T(),
                    key: r,
                    saltLength: i
                }).verify(t, e), !0;
            } catch (t) {
                if (!(t instanceof n)) throw t;
            }
            return !1;
        }
    };
}({}, function() {
    return this;
}());
//# sourceMappingURL=asmcrypto.js.map