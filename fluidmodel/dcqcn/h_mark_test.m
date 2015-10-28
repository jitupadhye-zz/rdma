function out= h_mark_test(t, q, kmin, kmax)
global pmax;
out = 0;
if q <= kmin
      out = 0;
end
if q > kmax
      out = 1;
end
if q>kmin && q <= kmax
      
      out = (q-kmin)/(kmax-kmin)*pmax;
      %{
      p = (q-kmin)/(kmax-kmin)*pmax;
      if rand()<p
          out = 1;
      else
          out = 0;
      end
      %}
end

if t<=0
    out = 1;
end
