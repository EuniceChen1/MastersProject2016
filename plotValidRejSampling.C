{
  static_cast<TH1D*>(_file0->Get("Projected"))->Draw("EHIST");
  static_cast<TH1D*>(_file0->Get("Thrown"))->Draw("EHIST SAME");
}
