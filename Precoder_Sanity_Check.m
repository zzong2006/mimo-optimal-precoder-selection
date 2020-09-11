for i = 1:100
   F0_val = large_label(i, 1) ;
   F2_val = large_label(i, 3) ;
   F8_val = large_label(i, 9) ;
   F10_val = large_label(i, 11);
   
   
   if F0_val ~= F2_val
       fprintf('%d th Data : F0 (%.20f) and F2 (%.20f) is different.\n', i, F0_val, F2_val);
   end
   if F0_val ~= F8_val
       fprintf('%d th Data : F0 (%.20f) and F9 (%.20f) is different.\n', i, F0_val, F8_val);
   end
   if F0_val ~= F10_val
       fprintf('%d th Data : F0 (%.20f) and F10 (%.20f) is different.\n', i, F0_val, F10_val);
   end
   
   F1_val = large_label(i, 2) ;
   F3_val = large_label(i, 4) ;
   F9_val = large_label(i, 10) ;
   F11_val = large_label(i, 12);
   
   if F1_val ~= F3_val
       fprintf('%d th Data : F1 (%.20f) and F3 (%.20f) is different.\n', i, F1_val, F3_val);
   end
   if F1_val ~= F9_val
       fprintf('%d th Data : F1 (%.20f) and F9 (%.20f) is different.\n', i, F1_val, F9_val);
   end
   if F1_val ~= F11_val
       fprintf('%d th Data : F1 (%.20f) and F11 (%.20f) is different.\n', i, F1_val, F11_val);
   end
   
   F4_val = large_label(i, 5);
   F6_val = large_label(i, 7);
   
   if F4_val ~= F6_val
       fprintf('%d th Data : F4 (%.20f) and F6 (%.20f) is different.\n', i, F4_val, F6_val);
   end
   
   F5_val = large_label(i, 6);
   F7_val = large_label(i, 8);
   
   if F5_val ~= F7_val
       fprintf('%d th Data : F5 (%.20f) and F7 (%.20f) is different.\n', i, F5_val, F7_val);
   end
   
   F12_val = large_label(i, 13);
   F13_val = large_label(i, 14);
   F14_val = large_label(i, 15);
   F15_val = large_label(i, 16);
   
   if F12_val ~= F13_val
       fprintf('%d th Data : F12 (%.20f) and F13 (%.20f) is different.\n', i, F12_val, F13_val);
   end
   if F12_val ~= F14_val
       fprintf('%d th Data : F12 (%.20f) and F14 (%.20f) is different.\n', i, F12_val, F14_val);
   end
   if F12_val ~= F15_val
       fprintf('%d th Data : F12 (%.20f) and F15 (%.20f) is different.\n', i, F12_val, F15_val);
   end

end