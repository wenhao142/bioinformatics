import "./globals.css";
import { ReactNode } from "react";

export const metadata = {
  title: "AD Locus Evidence Platform",
  description: "AD multi-omics scaffold",
};

export default function RootLayout({ children }: { children: ReactNode }) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
